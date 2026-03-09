# GPU acceleration

Rastair 2.1 and newer support running its @ML inference on @GPU:pl.
This can speed up the `call` subcommand significantly.

## Usage

Pass `--gpu` when running `rastair call`.

## Requirements

The machine running Rastair will need to have access to a GPU.

On Linux, @Vulkan drivers are needed.
On macOS, the system-provided Metal drivers are used.

```admonish warning
Rastair does not use Nvidia's CUDA!
If you are running other @ML workloads, it might still be necessary to ensure Vulkan drivers are installed.
```

For example, when using @SLURM, you can request a job with access to a @GPU
using `-p gpu_short --gres gpu:1`.
You might also need to run `module load Vulkan;` before running Rastair.

### Tested Hardware

- Built-in GPU on M4 MacBook Pro
- Dedicated RTX 8000 Quadro
