# ============================================================
# Modelo gráfico probabilístico sobre grafo dirigido enriquecido
# para generar un mapa probabilístico de trayectoria ciclónica
#
# Entrada:
#   data/grafo_aristas.csv
#   data/nodos.csv
#
# Salidas:
#   results/pgm_matriz_transicion_condicionada.csv
#   results/pgm_trayectorias_simuladas.csv
#   results/pgm_mapa_probabilistico.csv
#   results/pgm_mapa_probabilistico_no_cero.csv
#   results/pgm_mapa_probabilistico.png
# ============================================================

import re
import math
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# 1. Configuración general
# ============================================================

ARCHIVO_ARISTAS = Path("data") / "grafo_aristas.csv"
ARCHIVO_NODOS = Path("data") / "nodos.csv"

CARPETA_SALIDA = Path("results")
CARPETA_SALIDA.mkdir(exist_ok=True)

# Columnas del archivo de aristas
COL_ORIGEN = "node_origen"
COL_DESTINO = "node_destino"
COL_FRECUENCIA = "frecuencia_transicion"
COL_PROB = "probabilidad_transicion"

# Parámetros de simulación
N_SIMULACIONES = 10000
HORIZONTE_PASOS = 8

# Nodo inicial.
# Si se deja como None, el modelo elige el nodo inicial según frecuencia histórica de salida.
START_NODE = None
# Ejemplo:
# START_NODE = "HDB_20"

# Resolución espacial del mapa.
# 0.10 grados equivale aproximadamente a 11 km en latitud.
RESOLUCION_GRADOS = 0.10

# Ancho del corredor alrededor de cada trayectoria simulada.
# Una celda se considera alcanzada si está dentro de este radio.
CORREDOR_KM = 35.0

# Buffer para extender el dominio del mapa alrededor de los nodos.
BUFFER_GRADOS = 1.0

# Semilla para reproducibilidad
SEED = 42

EVIDENCIA = {
    # Condicionamiento temporal por mes.
    # Usa la columna "meses_donde_ocurre", por ejemplo: "8 (9), 9 (6), 7 (5)".
    # "mes": 9,

    # Ejemplo de condicionamiento por velocidad:
    # "velocidad": {
    #     "col": "velocidad_media",
    #     "value": 25.0,
    #     "sigma": 7.0
    # },

    # Ejemplo de condicionamiento por distancia a tierra:
    # "distancia_tierra": {
    #     "col": "cem_distancia_tierra_km_media_arista",
    #     "value": 250.0,
    #     "sigma": 150.0
    # },

    # Ejemplo de condicionamiento orográfico:
    # "orografia_250km": {
    #     "col": "cem_interaccion_orografica_250km_media_arista",
    #     "value": 0.25,
    #     "sigma": 0.20
    # }
}


# ============================================================
# 2. Funciones auxiliares
# ============================================================

def normalizar_vector(x):
    """
    Normaliza un vector para que sume 1.
    Si la suma es cero, regresa distribución uniforme.
    """
    x = np.asarray(x, dtype=float)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    x[x < 0] = 0

    suma = x.sum()

    if suma <= 0:
        return np.ones_like(x) / len(x)

    return x / suma


def haversine_km(lat1, lon1, lat2, lon2):
    """
    Distancia haversine en kilómetros.
    Acepta escalares o arreglos numpy.
    """
    R = 6371.0

    lat1 = np.radians(lat1)
    lon1 = np.radians(lon1)
    lat2 = np.radians(lat2)
    lon2 = np.radians(lon2)

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    )

    c = 2.0 * np.arcsin(np.sqrt(a))

    return R * c


def parsear_meses(texto):
    """
    Convierte una cadena como:
        "8 (9), 9 (6), 7 (5), 10 (2)"
    en:
        {8: 9, 9: 6, 7: 5, 10: 2}
    """
    if pd.isna(texto):
        return {}

    texto = str(texto)
    pares = re.findall(r"(\d+)\s*\((\d+)\)", texto)

    conteos = {}

    for mes, conteo in pares:
        mes = int(mes)
        conteo = int(conteo)
        conteos[mes] = conteos.get(mes, 0) + conteo

    return conteos


def score_mes(row, mes_objetivo, alpha=0.5):
    """
    Calcula un score probabilístico para una arista dado un mes objetivo.

    Usa suavizado:
        score = (conteo_mes + alpha) / (total + 12 * alpha)
    """
    if "meses_donde_ocurre" not in row.index:
        return 1.0

    conteos = parsear_meses(row["meses_donde_ocurre"])

    if len(conteos) == 0:
        return 1.0

    total = sum(conteos.values())

    return (conteos.get(mes_objetivo, 0) + alpha) / (total + 12.0 * alpha)


def score_gaussiano(valor, objetivo, sigma):
    """
    Score tipo verosimilitud gaussiana.
    Mientras más cerca esté valor de objetivo, mayor será el score.
    """
    if pd.isna(valor):
        return 1.0

    if sigma <= 0:
        raise ValueError("sigma debe ser mayor que cero.")

    z = (float(valor) - float(objetivo)) / float(sigma)

    return math.exp(-0.5 * z * z)


def calcular_score_evidencia(row, evidencia):
    """
    Calcula el factor de ajuste de una arista dada cierta evidencia.

    Modelo:
        P(Z_{t+1}=j | Z_t=i, X=e)
        proporcional a:
        P_hist(j | i) * score(e | arista i->j)
    """
    score = 1.0

    # Condicionamiento por mes
    if "mes" in evidencia and evidencia["mes"] is not None:
        score *= score_mes(row, evidencia["mes"])

    # Condicionamientos continuos
    for nombre, config in evidencia.items():
        if nombre == "mes":
            continue

        if not isinstance(config, dict):
            continue

        col = config.get("col", None)
        value = config.get("value", None)
        sigma = config.get("sigma", None)

        if col is None or value is None or sigma is None:
            continue

        if col not in row.index:
            continue

        score *= score_gaussiano(row[col], value, sigma)

    return score


def leer_aristas(path):
    """
    Lee y valida el archivo de aristas.
    """
    if not path.exists():
        raise FileNotFoundError(f"No se encontró el archivo: {path}")

    df = pd.read_csv(path)

    columnas_requeridas = [
        COL_ORIGEN,
        COL_DESTINO,
        COL_FRECUENCIA,
        COL_PROB
    ]

    for col in columnas_requeridas:
        if col not in df.columns:
            raise ValueError(f"No se encontró la columna requerida: {col}")

    df[COL_ORIGEN] = df[COL_ORIGEN].astype(str)
    df[COL_DESTINO] = df[COL_DESTINO].astype(str)
    df[COL_FRECUENCIA] = pd.to_numeric(df[COL_FRECUENCIA], errors="coerce").fillna(0)
    df[COL_PROB] = pd.to_numeric(df[COL_PROB], errors="coerce").fillna(0)

    return df


def leer_nodos(path, nodos_requeridos):
    """
    Lee la tabla de nodos.

    Si no existe, genera una plantilla en data/nodos.csv
    para que el usuario agregue latitud y longitud.
    """
    if not path.exists():
        plantilla = pd.DataFrame({
            "node_id": sorted(nodos_requeridos),
            "latitud": np.nan,
            "longitud": np.nan
        })

        path.parent.mkdir(exist_ok=True)
        plantilla.to_csv(path, index=False)

        raise FileNotFoundError(
            "\nNo se encontró data/nodos.csv.\n"
            "Se generó una plantilla en data/nodos.csv.\n"
            "Llena las columnas latitud y longitud para cada nodo HDB "
            "y vuelve a ejecutar el script.\n"
        )

    nodos = pd.read_csv(path)

    # Permitir algunos nombres alternativos
    rename = {}

    if "lat" in nodos.columns and "latitud" not in nodos.columns:
        rename["lat"] = "latitud"

    if "lon" in nodos.columns and "longitud" not in nodos.columns:
        rename["lon"] = "longitud"

    if "longitude" in nodos.columns and "longitud" not in nodos.columns:
        rename["longitude"] = "longitud"

    if "latitude" in nodos.columns and "latitud" not in nodos.columns:
        rename["latitude"] = "latitud"

    if "id" in nodos.columns and "node_id" not in nodos.columns:
        rename["id"] = "node_id"

    nodos = nodos.rename(columns=rename)

    columnas_requeridas = ["node_id", "latitud", "longitud"]

    for col in columnas_requeridas:
        if col not in nodos.columns:
            raise ValueError(
                f"El archivo data/nodos.csv debe contener la columna: {col}"
            )

    nodos["node_id"] = nodos["node_id"].astype(str)
    nodos["latitud"] = pd.to_numeric(nodos["latitud"], errors="coerce")
    nodos["longitud"] = pd.to_numeric(nodos["longitud"], errors="coerce")

    if nodos[["latitud", "longitud"]].isna().any().any():
        raise ValueError(
            "Hay nodos sin latitud o longitud. Completa data/nodos.csv."
        )

    nodos_en_archivo = set(nodos["node_id"])
    faltantes = sorted(set(nodos_requeridos) - nodos_en_archivo)

    if len(faltantes) > 0:
        raise ValueError(
            "Faltan coordenadas para estos nodos en data/nodos.csv:\n"
            + ", ".join(faltantes)
        )

    return nodos


# ============================================================
# 3. Construcción del PGM sobre el grafo
# ============================================================

def construir_transiciones_condicionadas(df_aristas, evidencia):
    """
    Construye la matriz de transición condicionada.

    Base:
        P_hist(Z_{t+1}=j | Z_t=i)

    Condicionada:
        P(Z_{t+1}=j | Z_t=i, X=e)
        ∝ P_hist(j | i) * score(e | arista i->j)
    """
    df = df_aristas.copy()

    scores = []

    for _, row in df.iterrows():
        scores.append(calcular_score_evidencia(row, evidencia))

    df["score_evidencia"] = scores

    df["peso_no_normalizado"] = df[COL_PROB] * df["score_evidencia"]

    probs_condicionadas = []

    for origen, grupo in df.groupby(COL_ORIGEN):
        pesos = grupo["peso_no_normalizado"].values

        # Si la evidencia anuló todas las probabilidades,
        # se regresa al modelo histórico base para ese nodo.
        if np.nansum(pesos) <= 0:
            pesos = grupo[COL_PROB].values

        p_norm = normalizar_vector(pesos)
        probs_condicionadas.extend(p_norm)

    # Cuidado: groupby conserva los índices originales.
    # Por eso asignamos por grupos respetando índice.
    df["probabilidad_condicionada"] = np.nan

    for origen, grupo in df.groupby(COL_ORIGEN):
        pesos = grupo["peso_no_normalizado"].values

        if np.nansum(pesos) <= 0:
            pesos = grupo[COL_PROB].values

        p_norm = normalizar_vector(pesos)

        df.loc[grupo.index, "probabilidad_condicionada"] = p_norm

    matriz = df.pivot_table(
        index=COL_ORIGEN,
        columns=COL_DESTINO,
        values="probabilidad_condicionada",
        aggfunc="sum",
        fill_value=0
    )

    return df, matriz


def construir_grafo(df_transiciones):
    """
    Construye el grafo dirigido enriquecido con la probabilidad condicionada.
    """
    G = nx.DiGraph()

    for _, row in df_transiciones.iterrows():
        origen = row[COL_ORIGEN]
        destino = row[COL_DESTINO]

        attrs = row.dropna().to_dict()

        G.add_edge(origen, destino, **attrs)

        # Peso probabilístico para interpretación directa
        G[origen][destino]["weight"] = float(row["probabilidad_condicionada"])

        # Costo para rutas más probables
        p = float(row["probabilidad_condicionada"])
        G[origen][destino]["costo_probabilistico_condicionado"] = (
            -np.log(p) if p > 0 else np.inf
        )

    return G


def construir_diccionario_salidas(df_transiciones):
    """
    Construye un diccionario de salidas para simular cadenas de Markov.
    """
    salidas = {}

    for origen, grupo in df_transiciones.groupby(COL_ORIGEN):
        destinos = grupo[COL_DESTINO].astype(str).tolist()
        probs = grupo["probabilidad_condicionada"].values.astype(float)
        probs = normalizar_vector(probs)

        salidas[origen] = {
            "destinos": destinos,
            "probs": probs
        }

    return salidas


def distribucion_inicial(df_transiciones, start_node=None):
    """
    Define P(Z_0).

    Si start_node se especifica, P(Z_0=start_node)=1.
    Si no, se usa frecuencia histórica de salida.
    """
    if start_node is not None:
        return [start_node], np.array([1.0])

    frec = (
        df_transiciones
        .groupby(COL_ORIGEN)[COL_FRECUENCIA]
        .sum()
        .sort_index()
    )

    nodos = frec.index.astype(str).tolist()
    probs = normalizar_vector(frec.values)

    return nodos, probs


# ============================================================
# 4. Simulación Monte Carlo de trayectorias
# ============================================================

def simular_trayectorias(
    salidas,
    nodos_iniciales,
    probs_iniciales,
    n_simulaciones,
    horizonte_pasos,
    seed=42
):
    """
    Simula trayectorias sobre el grafo.

    Cada trayectoria es una realización de:
        Z_0, Z_1, ..., Z_H
    """
    rng = np.random.default_rng(seed)

    trayectorias = []

    for sim in range(n_simulaciones):
        actual = rng.choice(nodos_iniciales, p=probs_iniciales)

        nodos = [actual]
        aristas = []
        prob_trayectoria = 1.0

        for _ in range(horizonte_pasos):
            if actual not in salidas:
                break

            destinos = salidas[actual]["destinos"]
            probs = salidas[actual]["probs"]

            idx = rng.choice(len(destinos), p=probs)
            siguiente = destinos[idx]
            p = probs[idx]

            aristas.append((actual, siguiente))
            nodos.append(siguiente)

            prob_trayectoria *= p
            actual = siguiente

        trayectorias.append({
            "id_simulacion": sim,
            "nodos": nodos,
            "aristas": aristas,
            "longitud_pasos": len(aristas),
            "probabilidad_trayectoria": prob_trayectoria
        })

    return trayectorias


def guardar_trayectorias(trayectorias, path):
    """
    Guarda las trayectorias simuladas en CSV.
    """
    filas = []

    for tr in trayectorias:
        filas.append({
            "id_simulacion": tr["id_simulacion"],
            "longitud_pasos": tr["longitud_pasos"],
            "probabilidad_trayectoria": tr["probabilidad_trayectoria"],
            "secuencia_nodos": " -> ".join(tr["nodos"]),
            "secuencia_aristas": " | ".join(
                [f"{u}->{v}" for u, v in tr["aristas"]]
            )
        })

    pd.DataFrame(filas).to_csv(path, index=False)


# ============================================================
# 5. Rasterización probabilística de las trayectorias
# ============================================================

def interpolar_segmento(lat1, lon1, lat2, lon2, resolucion_km=10.0):
    """
    Interpola puntos entre dos nodos para representar una arista/corredor.
    """
    dist = haversine_km(lat1, lon1, lat2, lon2)

    n = max(2, int(np.ceil(dist / resolucion_km)))

    lats = np.linspace(lat1, lat2, n)
    lons = np.linspace(lon1, lon2, n)

    return lats, lons


def marcar_celdas_en_corredor(
    lat_punto,
    lon_punto,
    lat_values,
    lon_values,
    radio_km
):
    """
    Devuelve las celdas cercanas a un punto dentro de un radio dado.
    """
    lat_radius_deg = radio_km / 111.0

    cos_lat = max(0.2, abs(np.cos(np.radians(lat_punto))))
    lon_radius_deg = radio_km / (111.0 * cos_lat)

    i0 = np.searchsorted(lat_values, lat_punto - lat_radius_deg, side="left")
    i1 = np.searchsorted(lat_values, lat_punto + lat_radius_deg, side="right")

    j0 = np.searchsorted(lon_values, lon_punto - lon_radius_deg, side="left")
    j1 = np.searchsorted(lon_values, lon_punto + lon_radius_deg, side="right")

    i0 = max(i0, 0)
    j0 = max(j0, 0)
    i1 = min(i1, len(lat_values))
    j1 = min(j1, len(lon_values))

    if i0 >= i1 or j0 >= j1:
        return []

    sub_lats = lat_values[i0:i1]
    sub_lons = lon_values[j0:j1]

    LAT, LON = np.meshgrid(sub_lats, sub_lons, indexing="ij")

    dist = haversine_km(lat_punto, lon_punto, LAT, LON)

    ii, jj = np.where(dist <= radio_km)

    celdas = [(i0 + i, j0 + j) for i, j in zip(ii, jj)]

    return celdas


def generar_mapa_probabilistico(
    trayectorias,
    nodos_df,
    resolucion_grados,
    corredor_km,
    buffer_grados,
    n_simulaciones
):
    """
    Convierte trayectorias simuladas en una superficie espacial de probabilidad.

    La probabilidad de una celda se calcula como:
        número de simulaciones que tocaron la celda / número total de simulaciones

    Por construcción, los valores quedan entre 0 y 1.
    """
    coords = {
        row["node_id"]: (row["latitud"], row["longitud"])
        for _, row in nodos_df.iterrows()
    }

    lat_min = nodos_df["latitud"].min() - buffer_grados
    lat_max = nodos_df["latitud"].max() + buffer_grados
    lon_min = nodos_df["longitud"].min() - buffer_grados
    lon_max = nodos_df["longitud"].max() + buffer_grados

    lat_values = np.arange(lat_min, lat_max + resolucion_grados, resolucion_grados)
    lon_values = np.arange(lon_min, lon_max + resolucion_grados, resolucion_grados)

    hits = np.zeros((len(lat_values), len(lon_values)), dtype=int)

    for tr in trayectorias:
        celdas_visitadas = set()

        # Si la trayectoria no avanzó, marcamos el nodo inicial.
        if len(tr["aristas"]) == 0 and len(tr["nodos"]) > 0:
            nodo = tr["nodos"][0]

            if nodo in coords:
                lat, lon = coords[nodo]

                celdas = marcar_celdas_en_corredor(
                    lat,
                    lon,
                    lat_values,
                    lon_values,
                    corredor_km
                )

                celdas_visitadas.update(celdas)

        # Marcamos los corredores de cada arista de la trayectoria.
        for u, v in tr["aristas"]:
            if u not in coords or v not in coords:
                continue

            lat1, lon1 = coords[u]
            lat2, lon2 = coords[v]

            lats_seg, lons_seg = interpolar_segmento(
                lat1,
                lon1,
                lat2,
                lon2,
                resolucion_km=max(5.0, corredor_km / 2.0)
            )

            for lat_p, lon_p in zip(lats_seg, lons_seg):
                celdas = marcar_celdas_en_corredor(
                    lat_p,
                    lon_p,
                    lat_values,
                    lon_values,
                    corredor_km
                )

                celdas_visitadas.update(celdas)

        # Cada simulación aporta máximo 1 a cada celda.
        for i, j in celdas_visitadas:
            hits[i, j] += 1

    prob = hits / float(n_simulaciones)

    filas = []

    for i, lat in enumerate(lat_values):
        for j, lon in enumerate(lon_values):
            filas.append({
                "latitud": lat,
                "longitud": lon,
                "probabilidad": prob[i, j]
            })

    mapa_df = pd.DataFrame(filas)

    return mapa_df, prob, lat_values, lon_values


def graficar_mapa_probabilistico(
    prob,
    lat_values,
    lon_values,
    nodos_df,
    path_salida
):
    """
    Genera una visualización básica del mapa probabilístico.
    """
    plt.figure(figsize=(11, 9))

    extent = [
        lon_values.min(),
        lon_values.max(),
        lat_values.min(),
        lat_values.max()
    ]

    plt.imshow(
        prob,
        origin="lower",
        extent=extent,
        aspect="auto"
    )

    plt.colorbar(label="Probabilidad estimada de tránsito")

    plt.scatter(
        nodos_df["longitud"],
        nodos_df["latitud"],
        s=35,
        edgecolors="black"
    )

    for _, row in nodos_df.iterrows():
        plt.text(
            row["longitud"],
            row["latitud"],
            row["node_id"],
            fontsize=8,
            ha="left",
            va="bottom"
        )

    plt.xlabel("Longitud")
    plt.ylabel("Latitud")
    plt.title(
        "Mapa probabilístico de trayectoria ciclónica\n"
        "PGM sobre grafo dirigido enriquecido + Monte Carlo"
    )

    plt.tight_layout()
    plt.savefig(path_salida, dpi=300, bbox_inches="tight")
    plt.show()


# ============================================================
# 6. Ejecución principal
# ============================================================

def main():
    print("Leyendo archivo de aristas...")
    df_aristas = leer_aristas(ARCHIVO_ARISTAS)

    nodos_requeridos = sorted(
        set(df_aristas[COL_ORIGEN]).union(set(df_aristas[COL_DESTINO]))
    )

    print(f"Número de nodos detectados: {len(nodos_requeridos)}")
    print(f"Número de aristas detectadas: {len(df_aristas)}")

    print("\nLeyendo archivo de nodos...")
    nodos_df = leer_nodos(ARCHIVO_NODOS, nodos_requeridos)

    print("\nConstruyendo transiciones condicionadas...")
    df_transiciones, matriz_cond = construir_transiciones_condicionadas(
        df_aristas,
        EVIDENCIA
    )

    matriz_cond.to_csv(
        CARPETA_SALIDA / "pgm_matriz_transicion_condicionada.csv"
    )

    df_transiciones.to_csv(
        CARPETA_SALIDA / "pgm_aristas_condicionadas.csv",
        index=False
    )

    print("Matriz de transición condicionada guardada.")

    print("\nConstruyendo grafo dirigido enriquecido...")
    G = construir_grafo(df_transiciones)

    print(f"Nodos del grafo: {G.number_of_nodes()}")
    print(f"Aristas del grafo: {G.number_of_edges()}")

    nx.write_graphml(
        G,
        CARPETA_SALIDA / "pgm_grafo_dirigido_enriquecido.graphml"
    )

    nx.write_gexf(
        G,
        CARPETA_SALIDA / "pgm_grafo_dirigido_enriquecido.gexf"
    )

    print("Grafo exportado en GraphML y GEXF.")

    print("\nPreparando simulador Monte Carlo...")
    salidas = construir_diccionario_salidas(df_transiciones)

    nodos_ini, probs_ini = distribucion_inicial(
        df_transiciones,
        start_node=START_NODE
    )

    print("Distribución inicial:")
    for nodo, p in zip(nodos_ini, probs_ini):
        print(f"  {nodo}: {p:.4f}")

    print("\nSimulando trayectorias...")
    trayectorias = simular_trayectorias(
        salidas=salidas,
        nodos_iniciales=nodos_ini,
        probs_iniciales=probs_ini,
        n_simulaciones=N_SIMULACIONES,
        horizonte_pasos=HORIZONTE_PASOS,
        seed=SEED
    )

    guardar_trayectorias(
        trayectorias,
        CARPETA_SALIDA / "pgm_trayectorias_simuladas.csv"
    )

    print("Trayectorias simuladas guardadas.")

    print("\nGenerando mapa probabilístico...")
    mapa_df, prob, lat_values, lon_values = generar_mapa_probabilistico(
        trayectorias=trayectorias,
        nodos_df=nodos_df,
        resolucion_grados=RESOLUCION_GRADOS,
        corredor_km=CORREDOR_KM,
        buffer_grados=BUFFER_GRADOS,
        n_simulaciones=N_SIMULACIONES
    )

    mapa_df.to_csv(
        CARPETA_SALIDA / "pgm_mapa_probabilistico.csv",
        index=False
    )

    mapa_df[mapa_df["probabilidad"] > 0].to_csv(
        CARPETA_SALIDA / "pgm_mapa_probabilistico_no_cero.csv",
        index=False
    )

    print("Mapa probabilístico guardado como CSV.")

    print("\nGenerando visualización...")
    graficar_mapa_probabilistico(
        prob=prob,
        lat_values=lat_values,
        lon_values=lon_values,
        nodos_df=nodos_df,
        path_salida=CARPETA_SALIDA / "pgm_mapa_probabilistico.png"
    )

    print("\nProceso terminado correctamente.")
    print("\nArchivos generados:")
    print("- results/pgm_matriz_transicion_condicionada.csv")
    print("- results/pgm_aristas_condicionadas.csv")
    print("- results/pgm_trayectorias_simuladas.csv")
    print("- results/pgm_mapa_probabilistico.csv")
    print("- results/pgm_mapa_probabilistico_no_cero.csv")
    print("- results/pgm_mapa_probabilistico.png")
    print("- results/pgm_grafo_dirigido_enriquecido.graphml")
    print("- results/pgm_grafo_dirigido_enriquecido.gexf")


if __name__ == "__main__":
    main()