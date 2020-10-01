# cr-comp
This is a repository for the code developed to produce the results in the paper titled **"Cognitive radio network with coordinated multipoint joint transmission"** [DOI:10.1002/dac.3310](http://onlinelibrary.wiley.com/doi/10.1002/dac.3310/abstract)

### Abstract
Cognitive radio (CR) is considered to be a promising technology for future wireless networks to make opportunistic utilization of the unused or underused licensed spectrum. Meanwhile, coordinated multipoint joint transmission (CoMP JT) is another promising technique to improve the performance of cellular networks. In this paper, we propose a CR system with the CoMP JT technique. We develop an analytical model of the received signal-to-noise ratio at a CR to determine the energy detection threshold and the minimum number of required samples for energy detection–based spectrum sensing in a CR network (CRN) with CoMP JT technique. The performance of energy detection–based spectrum sensing under the developed analytical model is evaluated by simulation and found to be reliable. We formulate an optimization problem for a CRN with the CoMP JT technique to configure the channel allocation and user scheduling for maximizing the minimum throughput of the users. The problem is found to be a complex mixed-integer linear programming. We solve the problem using an optimization tool for several CRN instances by limiting the number of slots in frames. Further, we propose a heuristic-based simple channel allocation and user scheduling algorithm to maximize the minimum throughput of the users in CRNs with the CoMP JT technique. The proposed algorithm is evaluated via simulation and found to be very efficient.

### Code
The paper is organized into two parts. In the first part interference analysis of the proposed CRN with CoMP JT is performed. The code required to simulate the system model and figures 4 and 5 can be found in the `interference_analysis` folder. The second part considers the channel allocation and scheduling scheme and the corresponding code can be found in the `channel_allocation_scheduling` folder. To simulate the figures 7 to 12, the `main.m` file needs to be run.

There is a short description at the beginning of each file depicting the function of the code.


