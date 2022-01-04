BUILD_TPYE=RelWithDebInfo

cmake . -B build \
    -DCMAKE_TOOLCHAIN_FILE=E:/vcpkg/scripts/buildsystems/vcpkg.cmake \
    -DVCPKG_TARGET_TRIPLET=x64-windows-static \
    -DCMAKE_MSVC_RUNTIME_LIBRARY="MultiThreaded$<$<CONFIG:Debug>:Debug>"
    
cmake --build build --config $BUILD_TPYE

echo "================= Baseline =================="
./baseline

echo "================= Optimized ================="
./build/$BUILD_TPYE/main.exe

echo "================== Finish ==================="