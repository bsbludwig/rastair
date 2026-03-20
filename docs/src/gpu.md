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

### Tested Hardware

- Built-in GPU on M4 MacBook Pro
- NVIDIA RTX 6000/8000 Quadro
- NVIDIA A100