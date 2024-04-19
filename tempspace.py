import pandas as pd

# 读取CSV文件
df = pd.read_csv('D:\models\ecModel_batch.csv')
df_slice = df.loc[7175:8143]
# 假设需要处理的列名为'column_name'
# 使用正则表达式提取每个条目左侧的系数
df_slice['extracted_coefficient'] = df_slice['Reaction'].str.extract(r'^(\d+\.\d+)')

# 将提取的系数保存到新的CSV文件中
df_slice['extracted_coefficient'].to_csv(r'C:\Users\Liao\Desktop\coefficients.csv', index=False)
