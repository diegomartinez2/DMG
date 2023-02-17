#! python
# copyright (c) 2022 C.Y. Wong, myByways.com simplified Stable Diffusion v0.1

import os, sys, time
import torch
import numpy as np
from omegaconf import OmegaConf
from PIL import Image
from einops import rearrange
from pytorch_lightning import seed_everything
from contextlib import nullcontext
from ldm.util import instantiate_from_config
from ldm.models.diffusion.plms import PLMSSampler
from ldm.models.diffusion.ddim import DDIMSampler
from transformers import logging

PROMPTS = [         # --prompt, one or more in an array
    'A high definition cartoon of humans fighting aliens on mars',
]
NEGATIVES = [       # negative prompt, one or more, default None (or an empty array)
    'rover'
]

HEIGHT = 512        # --H, default 512, beyond causes M1 to crawl
WIDTH = 512         # --W, default 512, beyond causes M1 to crawl
FACTOR = 8          # --f downsampling factor, default 8

FIXED = 0           # --fixed_code, 1 for repeatable results, default 0
SEED = 42           # --seed, default 42
NOISE = 0.0         # --ddim_eta, 0 deterministic, no noise - 1.0 random noise, ignored for PLMS (must be 0)
PLMS = 0            # --plms, default 1 on M1 for txt2img but ignored for img2img (must be DDIM)
ITERATIONS = 1      # --n_iter, default 1
SCALE = 7.5         # --scale, 5 further -> 15 closer to prompt, default 7.5
STEPS = 50          # --ddim_steps, practically little improvement >50 but takes longer, default 50

IMAGE = None        # --init-img, img2img initial latent seed, default None
STRENGTH = 0.75     # --strength 0 more -> 1 less like image, default 0.75

FOLDER = 'outputs'  # --outdir for images and history file below
HISTORY = 'history.txt'
CONFIG = 'configs/stable-diffusion/v1-inference.yaml'
CHECKPOINT = 'models/ldm/stable-diffusion-v1/model.ckpt'

def seed_pre():
    if not FIXED:
        seed_everything(SEED)

def seed_post(device):
    if FIXED:
        seed_everything(SEED)
        return torch.randn([1, 4, HEIGHT // FACTOR, WIDTH // FACTOR], device='cpu').to(torch.device(device.type))
    return None

def load_model(config, ckpt=CHECKPOINT):
    pl_sd = torch.load(ckpt, map_location='cpu')
    sd = pl_sd['state_dict']
    model = instantiate_from_config(config.model)
    model.load_state_dict(sd, strict=False)
    return model

def set_device(model):
    if torch.backends.mps.is_available():
        device = torch.device('mps')
        precision = nullcontext
    elif torch.cuda.is_available():
        device = torch.device('cuda')
        precision = torch.autocast
    else:
        device = torch.device('cpu')
        precision = torch.autocast
    model.to(device.type)
    model.eval()
    return device, precision

def load_image(image_file):
    image = Image.open(image_file).convert('RGB')
    w, h = image.size
    w, h = map(lambda x: x - x % 32, (w, h))
    image = image.resize((w, h), resample=Image.Resampling.LANCZOS)
    image = np.array(image).astype(np.float32) / 255.0
    image = image[None].transpose(0, 3, 1, 2)
    image = torch.from_numpy(image)
    return 2.*image - 1.0

def setup_sampler(model):
    global NOISE
    if IMAGE:
        image = load_image(IMAGE).to(model.device.type)
        init_latent = model.get_first_stage_encoding(model.encode_first_stage(image))
        sampler = DDIMSampler(model)
        sampler.make_schedule(ddim_num_steps=STEPS, ddim_eta=NOISE, verbose=False)
        t_enc = int(STRENGTH * STEPS)
        sampler.t_enc = t_enc
        sampler.z_enc = sampler.stochastic_encode(init_latent, torch.tensor([t_enc]).to(model.device.type))
    elif PLMS:
        sampler = PLMSSampler(model)
        NOISE = 0
    else:
        sampler = DDIMSampler(model)
    return sampler

def get_prompts():
    global NEGATIVES
    if NEGATIVES is None:
        NEGATIVES = [''] * len(PROMPTS)
    else:
        NEGATIVES.extend([''] * (len(PROMPTS)-len(NEGATIVES)))
    return zip(PROMPTS, NEGATIVES)

def generate_samples(model, sampler, prompt, negative, start):
    uncond = model.get_learned_conditioning(negative) if SCALE != 1.0 else None
    cond = model.get_learned_conditioning(prompt)
    if IMAGE:
        samples = sampler.decode(sampler.z_enc, cond, sampler.t_enc,
            unconditional_guidance_scale=SCALE, unconditional_conditioning=uncond)
    else:
        shape = [4, HEIGHT // FACTOR, WIDTH // FACTOR]
        samples, _ = sampler.sample(S=STEPS, conditioning=cond, batch_size=1,
            shape=shape, verbose=False, unconditional_guidance_scale=SCALE,
            unconditional_conditioning=uncond, eta=NOISE, x_T=start)
    x_samples = model.decode_first_stage(samples)
    x_samples = torch.clamp((x_samples + 1.0) / 2.0, min=0.0, max=1.0)
    return x_samples

def save_image(image):
    name = f'{time.strftime("%Y%m%d_%H%M%S")}.png'
    image = 255. * rearrange(image.cpu().numpy(), 'c h w -> h w c')
    img = Image.fromarray(image.astype(np.uint8))
    img.save(os.path.join(FOLDER, name))
    return name

def save_history(name, prompt, negative):
    with open(os.path.join(FOLDER, HISTORY), 'a') as history:
        history.write(f'{name} -> {"PLMS" if PLMS else "DDIM"}, Seed={SEED}{" fixed" if FIXED else ""}, Scale={SCALE}, Steps={STEPS}, Noise={NOISE}')
        if IMAGE:
            history.write(f', Image={IMAGE}, Strength={STRENGTH}')
        if len(negative):
            history.write(f'\n + {prompt}\n - {negative}\n')
        else:
            history.write(f'\n + {prompt}\n')

def main():
    print('*** Loading Stable Diffusion - myByways.com simple-sd version 0.1')
    tic1 = time.time()
    logging.set_verbosity_error()
    os.makedirs(FOLDER, exist_ok=True)

    seed_pre()
    config = OmegaConf.load(CONFIG)
    model = load_model(config)
    device, precision_scope = set_device(model)
    sampler = setup_sampler(model)
    start_code = seed_post(device)

    toc1 = time.time()
    print(f'*** Model setup time: {(toc1 - tic1):.2f}s')

    counter = 0
    with torch.no_grad():
        with precision_scope(device.type):
            with model.ema_scope():

                for iteration in range(ITERATIONS):
                    for prompt, negative in get_prompts():
                        print(f'*** Iteration {iteration + 1}: {prompt}')
                        tic2 = time.time()
                        images = generate_samples(model, sampler, prompt, negative, start_code)
                        for image in images:
                            name = save_image(image)
                            save_history(name, prompt, negative)
                            print(f'*** Saved image: {name}')
                        toc2 = time.time()

                        print(f'*** Synthesis time: {(toc2 - tic2):.2f}s')
                        counter += len(images)

    print(f'*** Total time: {(toc2 - tic1):.2f}s')
    print(f'*** Saved {counter} image(s) to {FOLDER} folder.')

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print('*** User abort, goodbye.')
    except FileNotFoundError as e:
        print(f'*** {e}')
