cmake_minimum_required(VERSION 3.18)

project(OpenAcoustics
  DESCRIPTION "Room-acoustics solver kernels for OpenFDTD-X (ADR-0007 contract): FDTD (low band) + geometric acoustics (high band)"
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

# ── 浮動小数点の縮約 (FMA) を禁止する ─────────────────────────────
# GCC は既定 (-std=gnu11) で文をまたいで a*b + c を FMA 1 命令に縮約する
# (-ffp-contract=fast)。FMA を持つ命令セット向けのビルドでは、数学的に同じ
# 式でもループごとに縮約の有無が変わり、最後の 1 bit が変わる。外壁と
# 障害物面のビット一致 (判定 o — 不変条件 9) と「同じ入力・同じバイナリなら
# 同じ結果」の主張はこれに依存するので縮約させない (Ubuntu 26.04 / GCC 15
# の CI ランナーで (o) が落ちたのを、手元の -march=x86-64-v3 ビルドで再現し
# この指定で一致することを確認した)。MSVC は /fp:precise (既定) で縮約しない。
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    add_compile_options(-ffp-contract=off)
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

# ── 幾何音響ソルバー (高域担当) ───────────────────────────────────
# FDTD (低域) とは完全に独立した別バイナリ。共有するのは src/wav.c だけで、
# そのために必要な ac_fopen は src/ga_main.c が定義する (main.c は使わない)。
# OpenMP は使わない : レイのエコーグラム加算はリダクションになり、
# スレッド数によるビット一致を壊すため (CLAUDE.md の決定性の規則)。
add_executable(ofdx_acoustic_ga
    src/ga_main.c
    src/ga_input.c
    src/ga_ofdx.c
    src/ga_setup.c
    src/ga_trace.c
    src/ga_shape.c
    src/ga_synth.c
    src/ga_json.c
    src/wav.c
)

target_include_directories(ofdx_acoustic_ga PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(ofdx_acoustic_ga PRIVATE ${MATH_LIB})
