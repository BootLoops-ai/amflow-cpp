# Docker

> 英文版本: [DOCKER.md](DOCKER.md)

`ghcr.io/chang18/amflow-cpp` 把 **FLINT 3.4 + Kira 3.1 + FireFly + amflow_cli** 打包在单个 Linux 镜像里, 这样你不用自己编译依赖栈. 本页讲怎么用.

> **本镜像*不*含:**
> - **Fermat** (闭源, 学术用户免费, 但 license 禁止重新分发). 运行时必须用 bind-mount 把本地 Fermat 安装挂进去 — 见下文 ["Fermat: 必须 bind-mount"](#fermat-必须-bind-mount).
> - **Mathematica** (闭源商业软件). AMFlow.cpp 运行时**不**需要 MMA; MMA 只在维护者生成参考 oracle 值时用.

---

## 快速开始

```bash
# 1. 拉镜像 (或 `docker build`, 见 "源码构建")
docker pull ghcr.io/chang18/amflow-cpp:latest

# 2. 跑一个样例 bench: 1-loop cutbubble at eps = 1/100
docker run --rm \
    -v "$HOME/Ferl7:/usr/share/Ferl7:ro" \
    -v "$(pwd):/work" \
    ghcr.io/chang18/amflow-cpp:latest \
    /opt/amflow-cpp/tools/bench/cutbubble_1L_eps001_black_box_amflow_cpp.json \
    /work/out.json

# 3. 查看结果
cat out.json | python3 -m json.tool | head -40
```

首次运行需要拉几百 MB 的镜像; 后续运行会复用 cache 层.

---

## 打包了什么

| 组件 | 版本 | 镜像内路径 | License |
|---|---|---|---|
| `amflow_cli` | 1.1.0 | `/usr/local/bin/amflow_cli` | MIT |
| FLINT | 3.4.0 | `/opt/flint/{bin,include,lib}` | LGPL-3 |
| Kira | 3.1 | `/opt/kira/{bin,lib,share}` + `/usr/local/bin/kira` (symlink) | GPL-3 |
| FireFly | `kira-2` 分支 | `/opt/kira/lib/libfirefly.so.*` | (FireFly project license) |
| Oracle benches | 来自 v1.1.0 | `/opt/amflow-cpp/tools/bench/` | (项目 license) |

镜像基于 `ubuntu:22.04`; 运行时 stage 不带 `-dev` 包但保留了二进制需要的每个共享库 (`libginac11`, `libcln6`, `libyaml-cpp0.7`, `libjemalloc2`, `libgmp10`, `libmpfr6`, `zlib1g`, `libstdc++6`).

> Kira 用 `-Dkyotocabinet=false` 构建. Kyoto Cabinet 是 pyRed 的可选键值 store 后端, AMFlow.cpp 不用; 此外 Kira 3.1 的 `keyvaluedb.h` 在 kyoto-cabinet 分支上有一个 include 顺序问题 (`#include <kcpolydb.h>` 放在 `#include "pyred/config.h"` 之前), 在守护宏还没可见时会让编译失败, 所以我们关掉它.

---

## Fermat: 必须 bind-mount

Robert Lewis 的 Fermat 学术用户免费, 但 license 不允许重新分发. AMFlow.cpp 的 IBP 后端 (Kira) 需要它做多变量 rational 代数, 所以镜像**必须**在运行时看到你本地的 Fermat 安装.

镜像内默认的 Fermat 查找路径是 `/usr/share/Ferl7/fer64` (这跟 [`src/api/run_json.cpp`](../src/api/run_json.cpp) 里的硬编码默认值对应). 两种挂载方式可用:

**方式 A — 匹配默认路径** (推荐, 不需要 env var):

```bash
docker run --rm \
    -v /your/local/Ferl7:/usr/share/Ferl7:ro \
    -v "$(pwd):/work" \
    ghcr.io/chang18/amflow-cpp:latest \
    input.json /work/out.json
```

**方式 B — 自定义挂载 + `FERMATPATH` env var:**

```bash
docker run --rm \
    -v /custom/path:/opt/fermat:ro \
    -e FERMATPATH=/opt/fermat/fer64 \
    -v "$(pwd):/work" \
    ghcr.io/chang18/amflow-cpp:latest \
    input.json /work/out.json
```

Fermat 缺失时, 入口脚本在 Kira 被调用之前就以 exit code 2 退出, 附带清晰错误信息.

可从 <https://home.bway.net/lewis/> 拿 Fermat; 把 `Ferl7/` 目录复制到宿主任何位置, 然后 bind-mount 过去.

---

## 挂载语义

- `-v /your/Ferl7:/usr/share/Ferl7:ro` — Fermat 安装 (只读 OK).
- `-v "$(pwd):/work"` — 你的工作目录. 容器的 `WORKDIR` 是 `/work`, 所以 bench JSON 里的**相对**路径都会解到这里. 你写到 `/work/...` 的输出立刻在宿主可见.
- 输入文件路径: 传**绝对路径** (挂载了就用 `/work/your-input.json`, 自带 oracle 用 `/opt/amflow-cpp/tools/bench/...`).

bench JSON 里的 `work_dir` 是**容器内**路径. 自带 oracle 用 `/tmp/...` (ephemeral). 要持久化 Kira 中间状态, 把 `work_dir` 指到 `/work` 下面.

---

## 跑打包的 oracle bench

镜像打包了 v1.1.0 的所有 oracle bench. 列表:

```bash
docker run --rm ghcr.io/chang18/amflow-cpp:latest \
    ls /opt/amflow-cpp/tools/bench
```

一个 4-loop banana 示例 (~6 分钟, 与 MMA 量级相当):

```bash
docker run --rm \
    -v "$HOME/Ferl7:/usr/share/Ferl7:ro" \
    -v "$(pwd):/work" \
    ghcr.io/chang18/amflow-cpp:latest \
    /opt/amflow-cpp/tools/bench/banana_4loop_eps001_black_box_amflow_cpp.json \
    /work/banana_4loop_out.json
```

要对照 committed MMA 参考, 在宿主上跑比较脚本 (或在容器里用 Python 3 跑):

```bash
python3 tools/bench/compare_sampled.py \
    tools/bench/banana_4loop_eps001_mma_reference.json \
    banana_4loop_out.json
```

---

## 源码构建

想钉到特定 FLINT/Kira 版本或改 `Dockerfile` 时, 自己构建镜像:

```bash
git clone https://github.com/chang18/amflow-cpp
cd amflow-cpp
docker build -f docker/Dockerfile -t amflow-cpp:dev .
```

构建顺序跑三个大 stage:

1. **FLINT 3.4 源码构建** (~5–10 min). Ubuntu 22.04 自带 FLINT 2.8.4 太老; 我们从 GitHub 拉 tarball 用 `--disable-static` 编译以减小镜像体积.
2. **Kira 3.1 + FireFly (subproject)** (~10–15 min). `meson setup` 通过 `subprojects/firefly.wrap` 自动从 `gitlab.com/firefly-library/firefly` 拉 FireFly (`kira-2` revision). 我们设 `-Dflint=true` 让 FireFly 共享 FLINT 3.4 安装, 而不是拉 2.8.4 fallback.
3. **amflow_cli** (~2 min). 标准 `cmake + ninja`.

冷构建总 wallclock 因硬件而异, 约 20–30 分钟.

---

## 排错

### `ERROR: Fermat executable not found at "..."`

入口脚本找不到 Fermat. 仔细检查 `-v` 挂载 (见 ["Fermat: 必须 bind-mount"](#fermat-必须-bind-mount)). 容器内, `${FERMATPATH:-/usr/share/Ferl7/fer64}` 这个文件必须存在且有可执行权限.

### `kira::ExceptionInternal: ...` 或 Kira 崩溃

最常见原因是 **Fermat 版本不匹配** — 镜像基于某个 FireFly+Kira 组合构建, 该组合期望特定 Fermat 协议. 如果你本地 Fermat 是别的版本, 先试 `Ferl7` (项目测过的 baseline) 再报问题.

### 输出 JSON 空 / 输出文件 `permission denied`

Docker 默认以 `root` 写容器输出, 在宿主 bind mount 上可能写权限不够. 两种解决:

- **以你自己的用户跑:** 在 `docker run` flag 加 `--user "$(id -u):$(id -g)"`. 镜像不大, root vs user 没啥差别; 这纯粹是结果文件的归属问题.
- **让 root 拥有输出, 之后 `sudo chown`** 改回去.

### 拉镜像时 `unauthorized: authentication required`

镜像托管在 GitHub Container Registry (`ghcr.io`) 上, 是 public 的. 遇到 auth 错误通常是你的 docker client token 过期 — `docker logout ghcr.io` 一般能解决.

---

## 限制

- **目前仅 Linux/amd64.** arm64 需要 cross-compile (FLINT 和 Kira 在 ARM 上都能干净构建, 但 CI 还没接) 或原生 arm runner.
- **容器内不带 Mathematica 集成.** 想重新生成 `*_mma_reference.json` 对照上游 MMA 的话, 在宿主用你自己的 Mathematica 安装跑那些步骤 — C++ 端不需要 MMA, Docker 镜像也有意不带.
- **镜像约 1 GB.** 大头是 Kira + FireFly 编译产物; FLINT 约 120 MB. 瘦身需要要么更深的多 stage prune 拔掉没用的 FireFly object, 要么用 Distroless 基础镜像.

---

## 另见

- [`docs/USER_GUIDE.md`](USER_GUIDE.md) — JSON 输入格式, oracle 目录, 完整 CLI 选项.
- [`docs/FAQ.md`](FAQ.md) — 非 Docker 的 build / install 问题.
- [`tools/bench/README.md`](../tools/bench/README.md) — oracle bench 布局与再生流程.
