# plot_trends.py — УСТОЙЧИВАЯ ВЕРСИЯ (чистим CSV и строим тренды)
import os, glob, re, csv
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTDIR = "plots"
os.makedirs(OUTDIR, exist_ok=True)

# --- Поиск CSV ---
cands = []
for pat in ("results*.csv", "build/results*.csv", "**/results*.csv"):
    cands += glob.glob(pat, recursive=True)
if not cands:
    raise FileNotFoundError("Не найден results*.csv. Сначала запусти ./build/app")
INPUT = max(cands, key=lambda p: Path(p).stat().st_mtime)
print(f"Использую файл: {INPUT}")

# --- Пробуем угадать разделитель ---
raw_lines = Path(INPUT).read_text(encoding="utf-8", errors="ignore").splitlines()
sample = "\n".join(raw_lines[:10])
try:
    dialect = csv.Sniffer().sniff(sample, delimiters=";, \t")
    sep = dialect.delimiter
except Exception:
    sep = ";"

print(f"sep='{sep}' (десятичные определим после очистки)")

# --- Читаем КАК ТЕКСТ, чтобы потом очистить ---
df = pd.read_csv(INPUT, sep=sep, dtype=str, engine="python")

# Нормализуем заголовки
def normalize_col(c: str) -> str:
    c = (c or "").strip().lower().replace(" ", "")
    c = c.replace(",", "")  # "i,j,k" -> "ijk"
    return c

df.columns = [normalize_col(c) for c in df.columns]

# Приводим имя размера к 'n'
if 'n' not in df.columns:
    for cand in ('size','dim','nn'):
        if cand in df.columns:
            df = df.rename(columns={cand:'n'})
            break
if 'n' not in df.columns:
    # иногда пишут "matrixsize" и т.п.
    for c in list(df.columns):
        if c.startswith("n") or "size" in c:
            df = df.rename(columns={c: 'n'})
            break
if 'n' not in df.columns:
    raise ValueError(f"Не нашёл колонку размера среди: {list(df.columns)}")

# --- Очистка числовых ячеек ---
# заменим запятую на точку в числах и уберём лишние символы (например, 's')
num_keep = re.compile(r"[^0-9eE\+\-\.]")  # всё, что НЕ цифры/знак/точка/экспонента

for col in df.columns:
    # чистим каждую ячейку
    df[col] = (
        df[col]
        .astype(str)
        .str.replace(",", ".", regex=False)
        .str.replace(num_keep, "", regex=True)
        .replace({"": None})
    )

# Приведение к числам (ошибки -> NaN)
df['n'] = pd.to_numeric(df['n'], errors='coerce')
orders = [c for c in df.columns if c != 'n']
for c in orders:
    df[c] = pd.to_numeric(df[c], errors='coerce')

# Удалим строки без данных
df = df.dropna(subset=['n']).reset_index(drop=True)

# Если остались строки, где все измерения NaN — уберём
valid_mask = df[orders].notna().any(axis=1)
df = df[valid_mask].reset_index(drop=True)

# Диагностика
print("Колонки:", ['n'] + orders)
print("Строк после очистки:", len(df))
if len(df) == 0:
    raise RuntimeError("После очистки не осталось данных. Покажи первые 10 строк CSV — подстрою парсер.")

x = df['n'].to_numpy(dtype=float)

def r2_score(y, yhat):
    ss_res = np.sum((y - yhat)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    return 1.0 - ss_res/ss_tot if ss_tot > 0 else float('nan')

def fit_and_plot(order_name: str, y: np.ndarray):
    y = y.astype(float)
    # отбрасываем точки, где нет значения по этому порядку
    mask = ~np.isnan(y)
    xs, ys = x[mask], y[mask]
    if len(xs) < 2:
        print(f"[WARN] Недостаточно точек для {order_name} (len={len(xs)}). Пропускаю.")
        return None, None, None

    deg = min(3, max(1, len(xs)-1))  # чтобы polyfit не падал на 2-3 точках
    coeffs = np.polyfit(xs, ys, deg=deg)
    p = np.poly1d(coeffs)
    yhat = p(xs)
    r2 = r2_score(ys, yhat)

    x_min, x_max = xs.min(), xs.max()
    xs_smooth = np.linspace(x_min, x_max, 400)
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

report = []
for col in orders:
    res = fit_and_plot(col, df[col].to_numpy())
    if res[0] is not None:
        report.append((col, *res))

# Общий средний тренд (по доступным порядкам)
if report:
    present_cols = [name for (name, *_rest) in report]
    yavg = df[present_cols].astype(float).mean(axis=1).to_numpy()
    res = fit_and_plot("avg", yavg)
    if res[0] is not None:
        report.append(("avg", *res))

md_path = os.path.join(OUTDIR, "trend_equations.md")
with open(md_path, "w", encoding="utf-8") as f:
    f.write("# Тренды (аппроксимация)\n\n")
    for name, coeffs, r2, path in report:
        f.write(f"**{name}**:\n\n")
        if len(coeffs) == 4:
            f.write(f"- y = {coeffs[0]:.6e}·n³ + {coeffs[1]:.6e}·n² + {coeffs[2]:.6e}·n + {coeffs[3]:.6e}\n")
        elif len(coeffs) == 3:
            f.write(f"- y = {coeffs[0]:.6e}·n² + {coeffs[1]:.6e}·n + {coeffs[2]:.6e}\n")
        else:
            f.write(f"- y = {coeffs[0]:.6e}·n + {coeffs[1]:.6e}\n")
        f.write(f"- R² = {r2:.6f}\n")
        f.write(f"- График: {path}\n\n")

print(f"Сохранено {len(report)} графиков → {OUTDIR}/")
print(f"Markdown с формулами: {md_path}")
