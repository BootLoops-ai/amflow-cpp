# 参与贡献

> 英文版本: [CONTRIBUTING.md](CONTRIBUTING.md)

AMFlow.cpp 是一个由 parity 驱动的 AMFlow 算法 C++17 重实现. 做行为性变更之前, 请先读 [`docs/CONTRIBUTING.md`](docs/CONTRIBUTING.md) 里项目专用的工作流.

简要版:

- 上游 Mathematica AMFlow (<https://gitlab.com/multiloop-pku/amflow>) 是行为规范.
- 凡是用来确立或保护 parity 的 Mathematica 参考输出与原始 cache, 都要保留.
- 给改动的 domain 加聚焦测试 (`tests/test_<domain>_*.cpp`).
- 项目状态变化时同步更新 `README.md` 和 `docs/AUDIT_MMA_PARITY.md`.
- 不要靠放宽数值容差或拿掉 parity 检查来让 C++ 不一致 "通过".

本地验证:

```bash
cmake -S . -B build -DAMFLOW_BUILD_DRIVER=ON
cmake --build build -j
ctest --test-dir build --output-on-failure -j 4
```

依赖 Kira / Fermat 的测试在这些工具缺失时会 skip. Mathematica / Wolfram Engine 只在重新生成参考数据时需要.

## 关于本实现

本仓库的 C++17 源码主要由 AI 编程 agent 在维护者指导下开发. 每一项行为性变更都以"与上游 Mathematica 参考实现的数值一致性"作为闸门; `tools/bench/` 下的 228 oracle benchmark 与 548-case GoogleTest 测试套是任何改动都必须满足的契约.

## 联系方式

bug 反馈与功能请求请到 GitHub issue: <https://github.com/chang18/amflow-cpp/issues>.

合作咨询, 引用问题, 或其它不适合公开 issue 的事项, 请联系维护者 <3250800970@qq.com>.
