cmake_minimum_required(VERSION 3.18)

project(OpenAcoustics
  DESCRIPTION "Room-acoustics FDTD solver kernel for OpenFDTD-X (ADR-0007 contract)"
  LANGUAGES C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)   # 既定 Release (OpenPEEC と同じ)
endif()

# ── Windows (MSVC) 対応 (OpenPEEC portability.md を踏襲) ─────────────
# libm は存在しない (MATH_LIB を空にする)。ソースは UTF-8 (日本語コメント)。
# M_PI 等は _USE_MATH_DEFINES が必要 (本実装は自前の AC_PI を使う)。
if(MSVC)
    add_compile_options(/utf-8)
    add_compile_definitions(_USE_MATH_DEFINES _CRT_SECURE_NO_WARNINGS NOMINMAX)
    set(MATH_LIB "")
else()
    set(MATH_LIB m)
endif()

# OpenMP は任意 (無ければシリアルでビルド : macOS の AppleClang 等)。
# 並列はリダクションを持たないセル/面ループのみで、スレッド数によらず
# 結果はビット単位で一致する (acoustic_check.sh が判定する)。
option(WITH_OPENMP "Enable OpenMP parallel loops" ON)
if(WITH_OPENMP)
    find_package(OpenMP QUIET)
endif()

# output folder
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(ofdx_acoustic_fdtd
    src/main.c
    src/input_ofd.c
    src/input_ofdx.c
    src/fdtd.c
    src/wav.c
    src/jsonout.c
)

target_include_directories(ofdx_acoustic_fdtd PRIVATE ${PROJECT_SOURCE_DIR}/include)

if(WITH_OPENMP AND OpenMP_C_FOUND)
    target_link_libraries(ofdx_acoustic_fdtd PRIVATE OpenMP::OpenMP_C)
    message(STATUS "OpenMP が見つかりました。並列ビルドを有効にします。")
else()
    message(STATUS "OpenMP なし。シリアルでビルドします。")
endif()

target_link_libraries(ofdx_acoustic_fdtd PRIVATE ${MATH_LIB})
