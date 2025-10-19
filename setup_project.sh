# === УСТАНОВКА ВСЕГО ДЛЯ ПРОЕКТА test ===
echo "🔍 Проверка и установка зависимостей C++ и Python..."

# Обновление индексов пакетов
sudo apt update -y

# --- Проверка и установка C++ инструментов ---
for pkg in build-essential cmake g++; do
  if ! dpkg -s $pkg &>/dev/null; then
    echo "⬇️  Устанавливаю $pkg..."
    sudo apt install -y $pkg
  else
    echo "✅ $pkg уже установлен"
  fi
done

# --- Проверка Python и pip ---
if ! command -v python3 &>/dev/null; then
  echo "⬇️  Устанавливаю Python3..."
  sudo apt install -y python3
else
  echo "✅ Python3 уже установлен ($(python3 --version))"
fi

if ! command -v pip &>/dev/null; then
  echo "⬇️  Устанавливаю pip..."
  sudo apt install -y python3-pip
else
  echo "✅ pip уже установлен ($(pip --version))"
fi

# --- Проверка виртуального окружения ---
if [ ! -d ".venv" ]; then
  echo "⚙️  Создаю виртуальное окружение..."
  python3 -m venv .venv
else
  echo "✅ .venv уже существует"
fi

# --- Установка Python-библиотек ---
echo "📦 Проверяю/устанавливаю numpy, pandas, matplotlib..."
source .venv/bin/activate
pip install --upgrade pip >/dev/null
pip install numpy pandas matplotlib -q
deactivate

echo "✅ Всё готово!"
echo "Теперь можно запустить сборку проекта:"
echo "  cmake -B build -S . && cmake --build build -j && ./build/app"
echo "А затем анализ:"
echo "  source .venv/bin/activate && python3 plot_trends.py"
