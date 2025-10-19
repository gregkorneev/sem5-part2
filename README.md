sem5-part2

Проверка на наличие всех расширений: ./setup_project.sh

Сборка и запуск: rm -rf build && cmake -B build && cmake --build build && cd ./build && ./app && cd .. && python3 plot_trends.py
