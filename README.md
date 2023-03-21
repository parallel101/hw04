# 高性能并行编程与优化 - 第04讲的回家作业

通过 pull request 提交作业。会批分数，但是：

没有结业证书，回家作业仅仅作为评估学习效果和巩固知识的手段，不必为分数感到紧张 :)
量力而行，只要能在本课中，学到昨天的自己不懂的知识，就是胜利，没必要和别人攀比。
注意不要偷看别人的作业哦！

- 课件：https://github.com/parallel101/course
- 录播：https://space.bilibili.com/263032155

作业提交时间不限 :) 即使完结了还想交的话我也会看的~ 不过最好在下一讲开播前完成。

- 如何开 pull request：https://zhuanlan.zhihu.com/p/51199833
- 如何设置 https 代理：https://www.jianshu.com/p/b481d2a42274

## 评分规则

- 在你的电脑上加速了多少倍，就是多少分！请在 PR 描述中写明加速前后的用时数据。
- 最好详细解释一下为什么这样可以优化。会额外以乘法的形式加分。
- 比如你优化后加速了 50 倍，讲的很详细，所以分数乘 2，变成 100 分！
- 比如你优化后加速了 1000 倍，但是你的 PR 描述是空，所以分数乘 0，变成 0 分！

## 作业要求

利用这次课上所学知识，修改 main.cpp，优化其中的多体引力求解器：

- 不允许使用多线程并行
- 不允许做算法复杂度优化
- 可以针对编译器和平台优化，这次不要求跨平台
- 可以用 xmmintrin.h，如果你觉得编译器靠不住的话

初始数据:

Initial energy: -8.571526
Final energy: -8.511777
Time elapsed: 6646 ms

编译指令加入O3优化：

Initial energy: -8.571526
Final energy: -8.511777
Time elapsed: 1737 ms

将结构体OOP改成DOP

Initial energy: -8.571526
Final energy: -8.511777
Time elapsed: 1734 ms

加入编译指令

```
#pragma GCC ivdep
#pragma GCC unroll 4
```

Initial energy: -8.571302
Final energy: -8.511518
Time elapsed: 1587 ms

加上暴力火车头：

Initial energy: -8.571527
Final energy: -8.511723
Time elapsed: 1175 ms


加入编译指令：

`-ffast-math -march=native`

Initial energy: -8.571527
Final energy: -8.511747
Time elapsed: 210 ms
