def data_process():
    map = {}
    with open('../resource/MGI_DO.rpt.txt', 'r') as f:
        for line in f.readlines():
            line = line.strip()
            words = line.split('\t')
            if 'MGI:' in words[-1] and True if words[-2].isdigit() else False:
                map[words[-2] + '\t'] = words[-1] + '\t'

    with open('../resource/MouseNetV2_entrez.txt', '+r') as f:
        t = f.read()
        for k, v in map.items():
            t = t.replace(k, v)
        # 读写偏移位置移到最开始处
        f.seek(0, 0)
        f.write(t)


def print_data():
    str = ''
    with open('../resource/MouseNetV2_entrez.txt', 'r') as f:
        for line in f.readlines():
            words = line.split('\t')
            if 'MGI:' in words[0] and 'MGI:' in words[1]:
                str = str + line
    print(str)


def data_clear():
    str = ''
    with open('../resource/MGI_fnl.txt', '+r') as f:
        for line in f.readlines():
            # if line.split('\t')[0][0].isdigit() or line.split('\t')[1][0].isdigit():
            #     continue
            words = line.split('\t')
            str = str + words[0][4:] + '\t' + words[1][4:] + '\t' + words[2]
    print(str)


data_clear()
