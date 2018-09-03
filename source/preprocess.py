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


# 基因之间的相似性，浮点数
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


mp_graph = mp_process()
mgi_graph = mgi_process()
mgi2mp_graph = mgi2mp_process()


# PPG, PGG, PPGG, PGPG


def PP(mp):
    return mp_graph[mp] if mp in mp_graph.keys() else False


# 目前只考虑了一个表型只能对应一个基因的情况
def PG(mp):
    for key, value in mgi2mp_graph.items():
        if mp in value:
            return key
        else:
            continue
    return False


def GG(mgi):
    for key in mgi_graph.keys():
        if mgi in key:
            print(key)
            print(mgi_graph.get(key))
            return key, mgi_graph.get(key)
        else:
            continue
    return False


def GP(mgi):
    return mgi2mp_graph[mgi] if mgi in mgi2mp_graph.keys() else False


def PPG(mp):
    if PP(mp):
        for mp_2 in PP(mp):
            print(PG(mp_2))


def PGG(mp):
    if PG(mp):
        print(GG(PG(mp)))


def PPGG(mp):
    if PP(mp):
        for mp_2 in PP(mp):
            PGG(mp_2)


def PGPG(mp):
    if PG(mp):
        if GP(PG(mp)):
            print(PG(GP(PG(mp))))
# 0000003:0005375
# PPG('0000015')
# 0000111-1856276-
# PGG('0000111')
