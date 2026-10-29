@echo off
call "D:\Program Files\Microsoft Visual Studio\18\Community\Common7\Tools\VsDevCmd.bat" -arch=x64
if errorlevel 1 exit /b 1
set PATH=C:\Qt\6.12.0\msvc2022_64\bin;%PATH%
for %%R in (27fa7cd 44126f0 5f16670 c1551b0 caf9139 847b6e1) do (
    echo === Building and benchmarking %%R ===
    cmake -S out\benchmark\history\src\%%R -B out\benchmark\history\build\%%R -G Ninja -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=cl.exe -DCMAKE_PREFIX_PATH=C:/Qt/6.12.0/msvc2022_64
    if errorlevel 1 exit /b 1
    cmake --build out\benchmark\history\build\%%R --config Release --target ray1
    if errorlevel 1 exit /b 1
    copy /Y out\benchmark\history\build\%%R\ray1.exe out\benchmark\history\exec\ray1.exe
    if errorlevel 1 exit /b 1
    windeployqt --release --no-translations --compiler-runtime out\benchmark\history\exec\ray1.exe
    if errorlevel 1 exit /b 1
    out\benchmark\history\exec\ray1.exe --benchmark --warmup 30 --frames 180 --output out\benchmark\history\results\%%R.csv
    if errorlevel 1 exit /b 1
    copy /Y out\benchmark\history\exec\ray1.exe out\benchmark\history\exec\ray1-%%R.exe
)
exit /b 0
