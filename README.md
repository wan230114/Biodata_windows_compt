区域的滑动窗口统计
---

（一个示例，可用于后续区域滑动统计的参考模板）

- 1M为窗口，100kb为步长，进行滑动统计数据的峰值与平均值。




---
方法1：  
每100kb作为窗口开始，每统计完一个窗口输出结果。
1. 先生成窗口范围 {(chr1, 1, 10*6):[[datas1], [datas2], [datas3]], ...}
2. 逐行读取
3. 判断是否到达100kb
- 是，继续生成窗口，加入统计
- 否，只加入统计

方法2：（当前采用）  
对每一个点计算上下窗口，统计数据情况，最后计算完一条染色体按染色体导出。

---
待做
- [ ] 后续可以按位置导出, 进一步降低内存消耗
- [ ] 生成窗口范围: 需要判断是否到达末端
