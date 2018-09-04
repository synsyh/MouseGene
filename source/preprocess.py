import re


# 表型之间的相似性
def mp_process():
    mp_gragh = {}
    with open('../resource/test_MP', 'r') as f:
        content = f.read().strip()
        terms = content.split('[Term]')
        for term in terms:
            reg = re.compile(r'(?<=id: MP:)\d+')
            mp_match = reg.search(term)
            if mp_match:
                reg = re.compile(r'(?<=is_a: MP:)\d+')
                it = reg.finditer(term)
                if it:
                    mp_sim_list = []
                    for mp_sim_match in it:
                        mp_sim_list.append(mp_sim_match.group())
                    mp_gragh[mp_match.group(0)] = mp_sim_list
                else:
                    continue
    return mp_gragh


# 基因之间的相似性，浮点数，这个字典的key是两个基因number组成的tuple，value是两个基因的相似度
def mgi_process():
    mgi_gragh = {}
    with open('../resource/test_MGI_sim', 'r') as f:
        for line in f.readlines():
            list = line.strip().split('\t')
            if int(list[0]) < int(list[1]):
                tuple = (list[0], list[1])
                mgi_gragh[tuple] = list[2]
            else:
                tuple = (list[1], list[0])
                mgi_gragh[tuple] = list[-1]
    return mgi_gragh


# 基因对应的表型
def mgi2mp_process():
    mgi2mp_graph = {}
    mp_list = []
    with open('../resource/test_MGI', 'r') as f:
        for line in f.readlines():
            mgi_reg = re.compile(r'(?<=MGI:)\d+')
            mp_reg = re.compile(r'(?<=MP:)\d+')
            mgi_match = mgi_reg.search(line)
            mp_match = mp_reg.search(line)
            if mgi_match.group() in mgi2mp_graph.keys():
                mgi2mp_graph[mgi_match.group()].append(mp_match.group())
            else:
                mgi2mp_graph[mgi_match.group()] = []
                mgi2mp_graph[mgi_match.group()].append(mp_match.group())
            # mgi2mp_graph[mgi_match.group()].append(mp_match.group())
    return mgi2mp_graph


# 运行三个预处理程序，得到三个字典
mp_graph = mp_process()
mgi_graph = mgi_process()
mgi2mp_graph = mgi2mp_process()


# PP运算，输入表型mp，返回对应表型或False
def PP(mp):
    return mp_graph[mp] if mp in mp_graph.keys() else False


# PG预算，输入表型mp，返回对应基因或False，目前只考虑了一个表型只能对应一个基因的情况
def PG(mp):
    for key, value in mgi2mp_graph.items():
        if mp in value:
            return key
        else:
            continue
    return False


# GG运算，输入基因mgi，返回相似基因和相似度或False
def GG(mgi):
    for key in mgi_graph.keys():
        if mgi in key:
            return key, mgi_graph.get(key)
        else:
            continue
    return False


# GP运算，输入基因mgi，返回基因对应的表型或False
def GP(mgi):
    return mgi2mp_graph[mgi] if mgi in mgi2mp_graph.keys() else False


# PPG运算，先进行PP运算，然后对返回的每个mp进行PG运算
def PPG(mp):
    if PP(mp):
        for mp_2 in PP(mp):
            print(PG(mp_2))
    else:
        print('no MP conforming')


# PGG运算，先进行PG运算，再进行GG运算
def PGG(mp):
    pg_mp = PG(mp)
    if pg_mp:
        print(GG(pg_mp))
    else:
        print('no Gene conforming')


# PPGG运算，先进行PP运算，再进行PGG运算
def PPGG(mp):
    pp_mp = PP(mp)
    if pp_mp:
        for mp_2 in pp_mp:
            PGG(mp_2)


# PGPG运算，进行PG、GP、PG运算，
def PGPG(mp):
    pg_mp = PG(mp)
    if pg_mp:
        gp_pg_mp = GP(pg_mp)
        if gp_pg_mp:
            for mp_2 in gp_pg_mp:
                pg_mp_2 = PG(mp_2)
                if pg_mp_2 == pg_mp:
                    continue
                else:
                    if mp_2 == mp:
                        continue
                    else:
                        print(mp + '->' + pg_mp + '->' + mp_2 + '->' + pg_mp_2)
            print('complete')
        else:
            print('gene have no conforming mp')
    else:
        print('mp have no conforming Gene')


# PPG正确示例：0000003:0005375-1857473 输出基因
PPG('0000003')
# PGG正确示例：0000111-1856276-26451 输出基因和相似度
PGG('0000111')
# PPGG正确示例：0000005:0001778-1857473-26452 输出基因和相似度
PPGG('0000005')
# PGPG正确示例：0001778-1857473:0001779-1856276
PGPG('0001778')
