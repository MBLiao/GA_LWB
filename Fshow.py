from ModelOperation import *


def fs_fun(individual, model0, target_table, GeneProteinMap, tartgetRxn, ref, metrxn):
    with model0 as model:
        # 遍历基因编码
        for i, gene in enumerate(individual):
            if gene != 0:
                # 获取目标基因/酶的名称和调整方式
                target_gene = target_table.iloc[i]["gene_name"]

                target_enzyme = GeneProteinMap[target_gene]
                # print(target_gene, target_enzyme)
                # 根据操作调整模型
                if gene == 1:
                    print('E' * 12)
                    print('E' * 12)
                    enzymeUsage = ref.fluxes[target_enzyme]
                    if enzymeUsage <= 0.000000001:
                        model.reactions.get_by_id(target_enzyme).lower_bound = 0.000000004
                    else:
                        model.reactions.get_by_id(target_enzyme).lower_bound = enzymeUsage * 4
                    print(model.reactions.get_by_id(target_enzyme).bounds)
                elif gene == 2:
                    print('D' * 12)
                    print('D' * 12)
                    enzymeUsage = ref.fluxes[target_enzyme]
                    model.reactions.get_by_id(target_enzyme).upper_bound = enzymeUsage * 0.5
                    print(model.reactions.get_by_id(target_enzyme).bounds)

                elif gene == 3:
                    print('O' * 12)
                    print('O' * 12)
                    model.reactions.get_by_id(target_enzyme).upper_bound = 0
                    print(model.reactions.get_by_id(target_enzyme).bounds)

        try:

            ft1 = time.time()
            solution_m = MPMA.moma(model=model, solution=metrxn, linear=False)
            ft2 = time.time()

            status = solution_m.status
            if status == 'optimal':
                mutantYield_m = solution_m.fluxes[tartgetRxn] / solution_m.fluxes['r_1714_REV']
                print("---------------------------------------------MOMA----------------------------------------------")
                print("TIME: ", ft2 - ft1)
                growth = solution_m.fluxes['r_2111']
                product = solution_m.fluxes[tartgetRxn]
                glucose = solution_m.fluxes['r_1714_REV']
                print(growth, product, glucose, mutantYield_m)
                print('-----------------------------------------------------------------------------------------------')

            else:
                mutantYield_m = 0.01
        except:
            mutantYield_m = 0.01


    return mutantYield_m















