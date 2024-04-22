run 100, 200, 400, 800, 1600 cells on CPU
run 100, 200, ... cells on the GPU


tp = Cell throughput = ncells * tfinal / runtime = cell-ms/second
rt = runtime = tfinal / runtime

what to expect:

As the number of cells on CPU increases
- throughput will increase to a plateau
- it might drop a bit after that

As the number of cells on GPU increases
- throughput will increase to a plateau
- it will increase much faster, then hit a hard limit
- the plateau will be much higher

tp_gpu >> tp_cpu

We want the tp_gpu and tp_cpu to both be maximised - becase tp = tp_gpu + tp_cpu

We need rt of cpu and gpu to be the same - that way one does not wait for the other
- we are "back filling" the cpu to extend capabilities of the gpu
- so aim for rt on cpu higher than on the gpu, to ensure gpu runtime is maximised

run cpu and gpu only tests recording
* ncells
* tfinal
* runtime

plot
* ncells vs tp
* ncells vs rt
* rt vs tp

the rt vs tp plot for cpu and gpu should cross over at some point (cp
