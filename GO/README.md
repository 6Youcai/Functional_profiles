# 背景
为了查找某个研究领域的相关信息，生物学家往往要花费大量的时间，更糟糕的是，不同的生物学数据库可能会使用不同的术语，
好比是一些方言一样，这让信息查找更加麻烦，尤其是使得机器查找无章可循。Gene Ontology就是为了解决这种问题而发起的一个项目。

# 结构
Gene Ontology中最基本的概念是term。GO里面的每一个entry都有一个唯一的数字标记，形如GO:nnnnnnn，还有一个term名，
比如"cell", "fibroblast growth factor receptor binding"，或者"signal transduction"。每个term都属于一个ontology，
总共有三个ontology，它们分别是分子功能（molecular function）,细胞组分（cellular component）和生物学过程（biological process）。

- molecular function:
> 基因产物个体的功能，可以描述为分子水平的活性（activity），如催化（catalytic）或结合（binding）活性 

- cellular component
> 细胞组件本体论，亚细胞结构、位置和大分子复合物，细胞的每个部分和细胞外环境

- biological process
> 分子功能的有序组合，达成更广的生物功能，由一个或多个分子功能有序组合而产生的系列事件。

一个基因product可能会出现在不止一个cellular component里面，也可能会在很多biological process里面起作用，并且在其中发挥不同的molecular function。

Ontology中的term有两种相互关系，它们分别是`is_a`关系和`part_of`关系。is_a表示简单的包含关系，比如 nuclear chromosome is_a chromosome；
part_of稍微复杂，如nucleus part_of cell，核肯定是细胞的一部分，但有的细胞没有核。

# GO 数据库
GO数据库[http://www.geneontology.org/]，收集的是对各种物种基因功能进行限定和描述的标准词汇（term），
是国际标准化的基因功能描述分类系统。根据基因产物的相关生物学过程、细胞组分以及分子功能三个大类分别给予定义，
而每一大类下又包含更多层级具体term，这些定义与具体物种无关。
