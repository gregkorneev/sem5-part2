cmake_minimum_required(VERSION 3.12)
project(matmul_bench CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Если есть include/, добавим
if (EXISTS "${CMAKE_SOURCE_DIR}/include")
  include_directories(${CMAKE_SOURCE_DIR}/include)
endif()

# Найдём main.cpp: сначала в src/, потом в корне
if (EXISTS "${CMAKE_SOURCE_DIR}/src/main.cpp")
  set(MAIN_SRC "${CMAKE_SOURCE_DIR}/src/main.cpp")
elseif (EXISTS "${CMAKE_SOURCE_DIR}/main.cpp")
  set(MAIN_SRC "${CMAKE_SOURCE_DIR}/main.cpp")
else()
  message(FATAL_ERROR "Не найден main.cpp ни в src/, ни в корне проекта")
endif()

# Собираем исполняемый файл 'app'
add_executable(app
  ${MAIN_SRC}
)

# Оптимизации
if (MSVC)
  target_compile_options(app PRIVATE /O2 /DNDEBUG)
else()
  target_compile_options(app PRIVATE -O3 -DNDEBUG)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-march=native" HAS_MARCH_NATIVE)
  if (HAS_MARCH_NATIVE)
    target_compile_options(app PRIVATE -march=native)
  endif()
endif()
