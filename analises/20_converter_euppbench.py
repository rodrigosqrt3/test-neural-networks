# Converte a EUPPBench (estações; Demaeyer et al. 2023) em tabelas por passo de previsão,
# no formato da nossa pipeline: uma linha por (estação, rodada), com a observação de T2m e
# a média e o dp do ensemble de cada variável.
#
# Tarefa "R2F" de Feik et al. (2024): treino nos REFORECASTS (11 membros, duas rodadas por
# semana, 1997–2016) e teste nas PREVISÕES de 2017–2018 (51 membros, diárias), 24/72/120 h.
#
# Uso:  py analises/20_converter_euppbench.py [paises]      (padrão: os 5; ex.: austria)
# Saída: dados/processados/euppbench/eupp_{ref,fc}_h{024,072,120}.csv
#   colunas: station_id, init, valid, step, obs, crps_ens (CRPS do ensemble bruto),
#            <var>_media, <var>_dp para cada variável, t2m_min, t2m_max, metadados da estação.
#   Temperaturas em °C. CRPS do ensemble bruto pela fórmula com estatísticas de ordem.
import sys, os
import numpy as np, pandas as pd, xarray as xr

BASE = r"dados/brutos/euppbench"
SAIDA = r"dados/processados/euppbench"
PAISES = sys.argv[1:] or ["austria", "belgium", "france", "germany", "netherlands"]
PASSOS = [24, 72, 120]
KELVIN = {"t2m", "mx2t6", "mn2t6", "stl1", "t850"}
GRUPOS = {  # arquivo -> {variável no zarr: nome na tabela}
    "surface": {v: v for v in ["t2m", "tcc", "tcw", "tcwv", "u10", "v10", "u100", "v100", "cape", "cin",
                               "sd", "stl1", "swvl1", "vis"]},
    "surface_postprocessed": {v: v for v in ["cp6", "mn2t6", "mx2t6", "p10fg6", "slhf6", "sshf6", "ssr6",
                                             "ssrd6", "str6", "strd6", "tp6"]},
    "pressure_500": {"z": "z500"},
    "pressure_700": {"q": "q700", "u": "u700", "v": "v700"},
    "pressure_850": {"t": "t850"},
}

def abrir(nome):
    return xr.open_zarr(os.path.join(BASE, nome + ".zarr"), consolidated=True)

def crps_ensemble(x, y):
    # x: (..., m) membros, y: (...); CRPS = E|X - y| - E|X - X'|/2, com E|X - X'| por estatísticas de ordem
    x = np.sort(x, axis=-1); m = x.shape[-1]
    t1 = np.nanmean(np.abs(x - y[..., None]), axis=-1)
    i = np.arange(1, m + 1)
    t2 = 2 * np.sum(x * (2 * i - m - 1), axis=-1) / m**2
    return t1 - t2 / 2

def tabela(pais, tipo):
    # tipo "reforecasts": dims (station_id, time, number, year, step); "forecasts": sem year
    obs_nome = f"stations_{tipo}_observations_surface_{pais}"
    obs = abrir(obs_nome)
    partes = {p: {} for p in PASSOS}
    meta = None
    for grupo, vars_ in GRUPOS.items():
        ds = abrir(f"stations_ensemble_{tipo}_{grupo}_{pais}")
        passos = [type(ds.step.values[0].item())(p) for p in PASSOS]
        ds = ds.sel(step=passos)
        if meta is None:
            meta = ds
        for v, nome in vars_.items():
            a = ds[v].squeeze([d for d in ds[v].dims if d in ("surface", "isobaricInhPa", "depthBelowLandLayer")]).load()
            a = a.transpose("station_id", "time", *(["year"] if tipo == "reforecasts" else []), "step", "number").values
            if nome in KELVIN:
                a = a - 273.15
            for k, p in enumerate(PASSOS):
                ak = a[..., k, :]
                partes[p][nome + "_media"] = np.nanmean(ak, axis=-1)
                partes[p][nome + "_dp"] = np.nanstd(ak, axis=-1, ddof=1)
                if nome == "t2m":
                    partes[p]["t2m_min"] = np.nanmin(ak, axis=-1)
                    partes[p]["t2m_max"] = np.nanmax(ak, axis=-1)
                    partes[p]["_membros"] = ak
    st = meta.station_id.values
    est = pd.DataFrame({"station_id": st, "lat": meta.station_latitude.values, "lon": meta.station_longitude.values,
                        "altitude": meta.station_altitude.values, "model_altitude": meta.model_altitude.values,
                        "model_orography": meta.model_orography.values, "pais": pais})
    saidas = {}
    for k, p in enumerate(PASSOS):
        o = obs.sel(step=type(obs.step.values[0].item())(p))
        if tipo == "reforecasts":
            y = o.t2m.transpose("station_id", "time", "year").values - 273.15
            vt = meta.valid_time.sel(step=type(meta.step.values[0].item())(p)).transpose("time", "year").values
            valid = np.broadcast_to(vt, (len(st),) + vt.shape)
        else:
            y = o.t2m.transpose("station_id", "time").values - 273.15
            vt = (meta.time.values + np.timedelta64(p, "h"))
            valid = np.broadcast_to(vt, (len(st),) + vt.shape)
        d = {nome: arr.reshape(-1) for nome, arr in partes[p].items() if nome != "_membros"}
        d["obs"] = y.reshape(-1)
        d["crps_ens"] = crps_ensemble(partes[p]["_membros"], y).reshape(-1)
        d["station_id"] = np.broadcast_to(st.reshape((-1,) + (1,) * (y.ndim - 1)), y.shape).reshape(-1)
        d["valid"] = valid.reshape(-1)
        df = pd.DataFrame(d)
        df["init"] = df["valid"] - pd.Timedelta(hours=p)
        df["step"] = p
        df = df[np.isfinite(df["obs"]) & np.isfinite(df["t2m_media"])]
        saidas[p] = df.merge(est, on="station_id")
    return saidas

os.makedirs(SAIDA, exist_ok=True)
for tipo, rot in [("reforecasts", "ref"), ("forecasts", "fc")]:
    acum = {p: [] for p in PASSOS}
    for pais in PAISES:
        for p, df in tabela(pais, tipo).items():
            acum[p].append(df)
        print(tipo, pais, "ok", flush=True)
    for p in PASSOS:
        df = pd.concat(acum[p], ignore_index=True)
        arq = os.path.join(SAIDA, f"eupp_{rot}_h{p:03d}.csv")
        df.to_csv(arq, index=False, float_format="%.6g")
        print(arq, len(df), "linhas;", df.station_id.nunique(), "estações; CRPS ensemble bruto",
              round(float(df.crps_ens.mean()), 4), flush=True)
