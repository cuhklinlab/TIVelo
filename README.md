## TIVelo: RNA Velocity estimation utilizing trajectory inference at the cell cluster level

![Workflow](docs/fig/workflow_1.pdf)



RNA velocity inference is a valuable tool for understanding cell development, differentiation, and disease progression. However, existing RNA velocity inference methods typically rely on explicit assumptions of ordinary differential equations (ODE), which prohibits them to capture complex transcriptome expression patterns. In this study, we introduce TIVelo, a novel RNA velocity estimation approach that first determines the velocity direction at the cell cluster level based on trajectory inference, before estimating velocity for individual cells. TIVelo calculates an orientation score to infer the direction at the cluster level without an explicit ODE assumption, which effectively captures complex transcriptional patterns, avoiding potential inconsistencies in velocity estimation for genes that do not follow the simple ODE assumption. We validated the effectiveness of TIVelo by its application to 16 real datasets and the comparison with five benchmarking methods.

## Usage && Installation
Please follow the [Tutorials](https://tivelo.readthedocs.io/en/latest/) for installation and Usage.



