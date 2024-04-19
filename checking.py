import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 设置显示最大行数
pd.set_option('display.max_rows', None)  # 或者使用具体数字例如 pd.set_option('display.max_rows', 1000)

# 设置显示最大列数
pd.set_option('display.max_columns', None)  # 或者使用具体数字例如 pd.set_option('display.max_columns', 50)

# 设置显示的宽度，以防止行内容被截断
pd.set_option('display.width', None)  # 根据你的显示器或终端调整宽度

# 为了确保显示所有内容不被省略，可以增加列宽设置
pd.set_option('display.max_colwidth', None)

GeneProteinMap = {}

df = pd.read_excel(r'D:\models\ecModel_batch.xls')
print(df.describe())
GeneProteinMap = {gene: enz for gene, enz in zip(df.iloc[7175:8143, 3], df.iloc[7175:8143, 0])}
print(GeneProteinMap)
print(len(GeneProteinMap))


# target_table = pd.read_csv(r"C:\Users\Liao\Desktop\ecFSEOF_0.05_0.6_1.05_2PE_result.csv")
# print(len(target_table))
# print(target_table.describe())
# ge = []
# for gene in target_table['gene_name']:
#     if gene in GeneProteinMap:
#         continue
#     else:
#         ge.append(gene)
# print(ge, 'is not in Map')

# io = r"C:\Users\Liao\Desktop\Result_本科毕设\FLUX.xlsx"
# data = pd.read_excel(io)

# print("-" * 40)
# print(data.shape)
# varname = data.keys()
# print("-" * 40)
# print(varname)
# # 查看数据类型
# datatype = data.dtypes
# print("-" * 40)
# print(datatype)
# checkdata = data.isnull().sum()
# print("-" * 40)
# print(checkdata)
# # 或者采用info直接查看数据信息
# # print(data_copy.info())
# print(data.describe())
# print("-" * 40)
# zero_count = data.eq(0).sum()
# print(zero_count)
# print("-" * 40)

# non_zero_count = data.ne(0).sum()
# non_zero_count = non_zero_count.drop('Abbreviation')


files = [r"C:\Users\Liao\Desktop\Result_本科毕设\FVA_result\1000FVA_0.9_2_1.1.csv",
         r"C:\Users\Liao\Desktop\Result_本科毕设\FVA_result\1000FVA_0.9_2.csv",
         r"C:\Users\Liao\Desktop\Result_本科毕设\FVA_result\1000FVA_nothing.csv"]

# 存储结果
zero_counts = []
valid_counts = []

# 读取每个文件并计算
for file in files:
    data = pd.read_csv(file)
    filtered_data = data[data['Abbreviation'].str.startswith('draw_')]
    print(data.describe())
    zero_count = (data['difference'] == 0).sum()
    zero_counts.append(zero_count)
    # 计算'difference'列中有效数据（不为0且不为"#VALUE!"）的数量
    valid_count = data['difference'].apply(lambda x: x != 0 and x != "#VALUE!").sum()
    valid_counts.append(valid_count)


colors_o = [(33, 26, 62), (53, 38, 96), (85, 59, 148), (152, 114, 202), (142, 141, 200), (165, 151, 182),
            (191, 217, 229), (246, 231, 237), (95, 32, 61)
            , (208, 108, 157)]
colors = tuple(tuple(c / 255.0 for c in color) for color in colors_o)

# 创建柱状图
# labels = ['FVA_0.9_1_1.1', 'FVA_0.9_2_1.1', 'FVA_0.9_2', 'FVA_0.1_0.1', 'FVA_nothing']
labels = ['FVA_0.9_2_1.1', 'FVA_0.9_2', 'FVA_nothing']

# 第一个图形：零值数量
plt.figure(figsize=(10, 8))
plt.bar(labels, zero_counts, color=colors, alpha=0.7)
plt.title('Zero Counts in FVA', size=18)
plt.xlabel('Conditions', size=16)
plt.ylabel('Zero Counts', size=16)
plt.xticks(rotation=15)  # 设置x轴标签的旋转角度
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()

# 第二个图形：有效值数量
plt.figure(figsize=(10, 8))
plt.bar(labels, valid_counts, color=colors, alpha=0.7)
plt.title('Valid Counts in FVA', size=18)
plt.xlabel('Conditions', size=16)
plt.ylabel('Valid Counts', size=16)
plt.xticks(rotation=15)  # 设置x轴标签的旋转角度
plt.grid(True, linestyle='--', alpha=0.5)
plt.show()






# 绘制柱状图
plt.figure(figsize=(10, 6))  # 设置图形的大小
non_zero_count.plot(kind='bar', color=colors)  # kind='bar' 创建柱状图，可以自定义颜色
plt.title('Non-Zero Counts in Each Method', size=18)  # 添加标题
plt.xlabel('Methods', size=16)  # 添加x轴标签
plt.ylabel('Zero Count', size=16)  # 添加y轴标签
plt.xticks(rotation=0)  # 设置x轴标签的旋转角度
plt.grid(True, linestyle='--', alpha=0.4)  # 添加网格线，可设置样式和透明度
plt.show()  # 显示图形















column_name_1 = 'yield_wild'
column_name_2 = 'yield_mutant_m'
column_name_3 = 'BPCY_m'
column_name_4 = 'distance'
column_name_5 = 'fold'

plt.figure(figsize=(10, 6))
position = data[column_name_1][1]
plt.axhline(y=position, color='red', linestyle='--')
plt.scatter(data.index, data[column_name_2], c='darkblue')
plt.xlabel('Genes(Enzyme)', fontsize=16)
plt.ylabel('Mutant Yield', fontsize=16)

# plt.title(f'Scatter Plot of {column_name}')
# plt.xlabel('N')
# plt.ylabel(column_name)
plt.show()

plt.figure(figsize=(10, 6))
plt.scatter(data.index, data[column_name_5], c='green')
plt.xlabel('Genes(Enzyme)', fontsize=16)
plt.ylabel('FOLD', fontsize=16)
plt.show()


plt.figure(figsize=(10, 6))
plt.scatter(data.index, data[column_name_4], c='darkblue')
plt.xlabel('Genes(Enzyme)', fontsize=16)
plt.ylabel('Distance', fontsize=16)
plt.show()


# 绘制箱形图
plt.figure(figsize=(6, 6))
boxprops = dict(linestyle='-', linewidth=2, color='red', facecolor='yellow')
sns.boxplot(y=data[column_name_2], boxprops=boxprops)
plt.title(f'Box Plot of {column_name_2}')
plt.show()

plt.figure(figsize=(6, 6))
boxprops = dict(linestyle='-', linewidth=2, color='red', facecolor='yellow')
sns.boxplot(y=data[column_name_3], boxprops=boxprops)
plt.title(f'Box Plot of {column_name_3}')
plt.show()









