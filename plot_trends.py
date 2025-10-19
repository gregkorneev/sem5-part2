# ===========================================
# Построение трендов по results*.csv
# - Автопоиск свежего CSV (включая build/)
# - Устойчивое чтение: чистим числа от 's', пробелов, запятых
# - Для каждого порядка: точки + полином (до степени 3) + R^2
# - Отчёт в plots/trend_equations.md
# ===========================================

import os, glob, re, csv
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # не открываем окна — сохраняем PNG
import matplotlib.pyplot as plt

OUTDIR = "plots"
os.makedirs(OUTDIR, exist_ok=True)

# --- Ищем самый свежий results*.csv ---
cands = []
for pat in ("results*.csv", "build/results*.csv", "**/results*.csv"):
    cands += glob.glob(pat, recursive=True)
if not cands:
    raise FileNotFoundError("Не найден results*.csv. Сначала запусти ./build/app")
INPUT = max(cands, key=lambda p: Path(p).stat().st_mtime)
print(f"Использую файл: {INPUT}")

# --- Пытаемся угадать разделитель (обычно ';') ---
raw_lines = Path(INPUT).read_text(encoding="utf-8", errors="ignore").splitlines()
sample = "\n".join(raw_lines[:10])
try:
    dialect = csv.Sniffer().sniff(sample, delimiters=";, \t")
    sep = dialect.delimiter
except Exception:
    sep = ";"
print(f"sep='{sep}' (десятичные обработаем вручную)")

# --- Читаем как строки, чтобы безопасно «почистить» значения ---
df = pd.read_csv(INPUT, sep=sep, dtype=str, engine="python")

# Нормализация заголовков: убираем пробелы, переводим в нижний регистр, "i,j,k" -> "ijk"
def normalize_col(c: str) -> str:
    c = (c or "").strip().lower().replace(" ", "")
    c = c.replace(",", "")
    return c

df.columns = [normalize_col(c) for c in df.columns]

# Приводим столбец размера к имени 'n'
if "n" not in df.columns:
    for cand in ("size", "dim", "nn"):
        if cand in df.columns:
            df = df.rename(columns={cand: "n"})
            break
if "n" not in df.columns:
    # последний шанс: если столбец выглядит как 'n...' или содержит 'size'
    for c in list(df.columns):
        if c.startswith("n") or "size" in c:
            df = df.rename(columns={c: "n"})
            break
if "n" not in df.columns:
    raise ValueError(f"Не нашёл колонку размера среди: {list(df.columns)}")

# --- Чистка чисел: меняем запятую на точку, удаляем всё, кроме 0-9 eE +-. ---
num_keep = re.compile(r"[^0-9eE\+\-\.]")

for col in df.columns:
    df[col] = (
        df[col]
        .astype(str)
        .str.replace(",", ".", regex=False)     # '0,123s' -> '0.123s'
        .str.replace(num_keep, "", regex=True)  # '0.123s' -> '0.123'
        .replace({"": None})
    )

# --- Приводим к числам ---
df["n"] = pd.to_numeric(df["n"], errors="coerce")
orders = [c for c in df.columns if c != "n"]
for c in orders:
    df[c] = pd.to_numeric(df[c], errors="coerce")

# Удаляем строки без n и строки, где нет ни одного измерения
df = df.dropna(subset=["n"]).reset_index(drop=True)
valid_mask = df[orders].notna().any(axis=1)
df = df[valid_mask].reset_index(drop=True)

print("Колонки:", ["n"] + orders)
print("Строк после очистки:", len(df))
if len(df) == 0:
    raise RuntimeError("После очистки не осталось данных.")

x = df["n"].to_numpy(dtype=float)

# --- Метрики и отрисовка ---
def r2_score(y, yhat):
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")

def fit_and_plot(order_name: str, y_raw: np.ndarray):
    # Оставляем только точки, где есть данных по этому порядку
    y = y_raw.astype(float)
    mask = ~np.isnan(y)
    xs, ys = x[mask], y[mask]
    if len(xs) < 2:
        print(f"[WARN] Недостаточно точек для {order_name} (len={len(xs)}). Пропуск.")
        return None, None, None

    # Степень полинома: 3 (кубика). Если точек мало, понижаем степень.
    deg = min(3, max(1, len(xs) - 1))
    coeffs = np.polyfit(xs, ys, deg=deg)
    p = np.poly1d(coeffs)
    yhat = p(xs)
    r2 = r2_score(ys, yhat)

    xs_smooth = np.linspace(xs.min(), xs.max(), 400)
    ys_smooth = p(xs_smooth)

    plt.figure()
    plt.scatter(xs, ys, label="Измерения")
    plt.plot(xs_smooth, ys_smooth, label=f"Аппроксимация (deg={deg})")
    if deg == 3:
        eq = f"y = {coeffs[0]:.3e}·n³ + {coeffs[1]:.3e}·n² + {coeffs[2]:.3e}·n + {coeffs[3]:.3e}"
    elif deg == 2:
        eq = f"y = {coeffs[0]:.3e}·n² + {coeffs[1]:.3e}·n + {coeffs[2]:.3e}"
    else:
        eq = f"y = {coeffs[0]:.3e}·n + {coeffs[1]:.3e}"
    plt.title(f"Порядок {order_name}\n{eq}\nR² = {r2:.4f}")
    plt.xlabel("n")
    plt.ylabel("t, сек")
    plt.legend()
    plt.tight_layout()
    out = os.path.join(OUTDIR, f"trend_{order_name}.png")
    plt.savefig(out, dpi=160)
    plt.close()
    return coeffs, r2, out

# Строим графики по всем порядкам
report = []
for col in orders:
    res = fit_and_plot(col, df[col].to_numpy())
    if res[0] is not None:
        report.append((col, *res))

# Общий «средний» тренд (avg)
if report:
    present_cols = [name for (name, *_rest) in report]
    yavg = df[present_cols].astype(float).mean(axis=1).to_numpy()
    res = fit_and_plot("avg", yavg)
    if res[0] is not None:
        report.append(("avg", *res))

# Markdown-отчёт: формулы + R^2 + файлы графиков
md_path = os.path.join(OUTDIR, "trend_equations.md")
with open(md_path, "w", encoding="utf-8") as f:
    f.write("# Тренды (аппроксимация)\n\n")
    for name, coeffs, r2, path in report:
        if len(coeffs) == 4:
            eq = f"y = {coeffs[0]:.6e}·n³ + {coeffs[1]:.6e}·n² + {coeffs[2]:.6e}·n + {coeffs[3]:.6e}"
        elif len(coeffs) == 3:
            eq = f"y = {coeffs[0]:.6e}·n² + {coeffs[1]:.6e}·n + {coeffs[2]:.6e}"
        else:
            eq = f"y = {coeffs[0]:.6e}·n + {coeffs[1]:.6e}"
        f.write(f"**{name}**\n\n- {eq}\n- R² = {r2:.6f}\n- График: {path}\n\n")

print(f"Сохранено {len(report)} графиков → {OUTDIR}/")
print(f"Markdown с формулами: {md_path}")
