import random

# DEFAULT PARAMETERS
PS = 200
IT = 50
CP = 0.95
MP = 0.05

class Instance:
    def __init__(self, n, state, constant, n_col):
        self.n = n
        self.state = state
        self.constant = constant
        self.n_row = len(state)
        self.n_col = n_col

    def __getitem__(self, i):
        return self.ap[i]

    def __len__(self):
        return len(self.ap)

def create_output(genes):
    output_file = open("output.txt", "w")
    data = genes[1]
    for t in data:
        out_data = str(t[0]) + " " + str(t[1]) + "\n"
        output_file.write(out_data)
    output_file.close()


def compute_score(gene, instance):
    matrix = [[0 for i in range(instance.n_col)] for j in range(instance.n_row)]
    num = 0
    for t in gene:
        for i in xrange(t[0], t[0] + instance.constant[num]["M"]):
            matrix[0][i] += 1
            i += 1
        for j in xrange(t[0] + instance.constant[num]["M"], t[1]):
            matrix[1][j] += 1
            j += 1
        for k in xrange(t[1], t[1] + instance.constant[num]["O"]):
            matrix[2][k] += 1
            k += 1
        num += 1
    score = 0
    # score_L = 0
    # score_G = 0
    # score_T = 0
    for i, j, k in zip(matrix[0], matrix[1], matrix[2]):
        if i > instance.state["L"]:
            # score_L = 0
            score += 1
        if j > instance.state["G"]:
            # score_G = 0
            score += 1
        if k > instance.state["T"]:
            # score_T = 0
            score += 1
    return score


def run(instance, ps = PS, pc = CP, mit = IT):
    def init_population():
        """Generate initial population from random shuffles of the tasks."""
        population = []
        for i in xrange(ps):
            gene = []
            for j in xrange(instance.n):
                gene_a = random.randint(0, instance.constant[j]["R"])
                gene_b = random.randint(gene_a + instance.constant[j]["M"] + instance.constant[j]["S"], gene_a + instance.constant[j]["M"] + instance.constant[j]["C"])
                gene.append([gene_a, gene_b])
            population.append(gene)
        # --------- for debug ---------
        # population.append([[0, 60],[10, 80],[50, 130],[70, 150]])
        # -----------------------------
        return population

    def crossover(p1, p2):
        nt = len(p1) # total number of tasks
        child = []
        for i in xrange(nt):
            if random.getrandbits(1):
                child.append(p1[i])
            else:
                child.append(p2[i])
        return child

    def mutation(p):
        nj = len(p)
        nt = len(p[0])
        i = random.randint(0, nj - 1)
        j = random.randint(0, nt - 1)
        mu = random.randint(-5, 5)
        select = random.random()
        # t_gate = p[i][1] - p[i][0]
        if select < 0.5:
            p[i][j] += mu
        # for i in xrange(nj):
        #     if select < 0.8:
        #         temp = p[i][0] + mu
        #         if 0 <= temp <=instance.constant[i]["R"]:
        #             p[i][0] += mu
        #         if p[i][0] + instance.constant[i]["M"] + instance.constant[i]["S"] <= p[i][1] <= p[i][0] + instance.constant[i]["M"] + instance.constant[i]["C"]:
        #             p[i][1] += mu
        #     else:
        #         p[i][0] = random.randint(0, instance.constant[i]["R"])
        #         p[i][1] = random.randint(p[i][0] + instance.constant[i]["M"] + instance.constant[i]["S"], p[i][0] + instance.constant[i]["M"] + instance.constant[i]["C"])

    def repair(p):
        nt = len(p)
        for i in xrange(nt):
            # check output A
            if not 0 <= p[i][0]:
                p[i][0] = 0
            if not p[i][0] <= instance.constant[i]["R"]:
                p[i][0] -= p[i][0] - instance.constant[i]["R"]
            # check output B
            if not p[i][0] + instance.constant[i]["M"] + instance.constant[i]["S"] <= p[i][1]:
                p[i][1] += p[i][0] + instance.constant[i]["M"] + instance.constant[i]["S"] - p[i][1]
            if not p[i][1] <= p[i][0] + instance.constant[i]["M"] + instance.constant[i]["C"]:
                p[i][1] -= p[i][1] - (p[i][0] + instance.constant[i]["M"] + instance.constant[i]["C"])
    while True:
        pm = MP
        child_pop = []
        pop = [(compute_score(pp, instance), pp) for pp in init_population()]
        pop.sort()
        print(pop)
        pre_score = pop[0][0]
        if pre_score == 0:
            return pop[0]
        local_max_counter = 0
        for it in xrange(mit):
            h_pop = len(pop) / 2
            # create the half number of genes for new population
            for n in xrange(h_pop):
                # Create two new elements for crossover
                elite_p = pop[random.randint(0, ps - 1)][1]
                rand_p = pop[random.randint(0, ps - 1)][1]
                if random.random() < pc:
                    child_p = crossover(elite_p, rand_p)
                    if random.random() < pm:
                        mutation(child_p)
                    repair(child_p)
                    child_pop.append((compute_score(child_p, instance), child_p))
            pop = pop + child_pop
            pop.sort()
            pop = pop[:ps]
            print(pop)
            if pop[0][0] == 0:
                return pop[0]
            # # check local_max
            if pop[0][0] == pre_score:
                local_max_counter += 1
                if local_max_counter == 10 and pre_score < 10:
                    pm = 0.8
            else:
                pre_score = pop[0][0]
                local_max_counter = 0
    return pop[0]


def create_instance(input):
    """
    create instance
    :param input: <list>
    :return: <Instance>
    """
    st_cap = [int(x) for x in input[0].rstrip().split(" ")]
    state = dict(zip(["L", "G", "T"], st_cap))
    key = ["R", "M", "S", "O", "C"]
    const = []

    max_time = 0
    for consts in input[2]:
        const_str = consts.rstrip().split(" ")
        const_int = [int (x) for x in const_str]
        const_dict = dict(zip(key,const_int))
        const.append(const_dict)
        max_time = max(max_time, sum(const_int))
    return Instance(int(input[1]), state, const, max_time)


def open_file():
    """
    open the text file and read data into a list
    :return: <list> - string
    """
    input_file = open("input3.txt")
    place = input_file.readline()
    num_airplanes = input_file.readline()
    airplanes = input_file.readlines()
    input_file.close()
    return place, num_airplanes, airplanes


def main():
    instance = create_instance(open_file())
    genes = run(instance, ps=PS, mit=IT, pc=CP)
    print genes
    create_output(genes)

# Program starts here
if __name__ == "__main__":
    main()