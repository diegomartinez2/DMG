#! python -O
# copyright (c) 2022 C.Y. Wong, myByways.com simplified Stable Diffusion v1.0
import os, time, pickle, warnings, logging, logging.handlers
import pickletools
import sys
os.environ['TRANSFORMERS_OFFLINE'] = '0' #transformers/utils/hub.py
import torch
import numpy as np
from pytorch_lightning import seed_everything
from omegaconf import OmegaConf
from PIL import Image

#######################################
# Global Variables
#######################################
import json

class CONFIG():
    """Global configuration variables - edit or override"""

    # Stable Diffusion `txt2img`:
    PROMPTS = []            # one or more prompts in an array, one per output image
    NEGATIVES = None        # negative prompt, one or more, default None / []
    TILING = False          # seamless tiling, default=False

    # Stable Diffusion `img2img` with mask and `inpaint`:
    IMAGE = None            # img2img initial seed image in PNG format, default None
    MASK = None             # either: inpainting mask in PNG format, white to overwrite, black to retain, default None / []
    MASK = []               # or: one or more prompts for CLIPSeg to generate the mask
    UNMASK = []             # and optionally: negative mask to deduct from MASK, default None / []

    # Real-ESRGAN upscaling and face enhancement:
    UPSCALE = 0             # default 0 to skip upscaling, or set to 1 for no size change but with enhancement
    FACE_ENHANCE = False    # when upscaling, enable face enhancement (runs on CPU only), default False

    # Also important:
    SEED = None             # Seed, default None for random (same seed for predictable txt2img results)
    WIDTH = 512             # image height default 512 (larger causes M1 to crawl)
    HEIGHT = 512            # image width default 512 (larger causes M1 to crawl)
    OUTPUT_DIR = 'outputs'          # output directory for images, generated masks, and history file
    HISTORY_FILE = 'history.txt'    # history log file name
    FILENAME = '%y%m%d_%H%M%S.png'

    class SD:               # Core Stable Diffusion:
        ITERATIONS = 1      # iterations (repeat runs using the same set of prompts), default 1
        FACTOR = 8          # down sampling factor, default 8 (can be larger, but not smaller)
        STEPS = 50          # number of steps, seldom visible improvement beyond ~50, default 50
        NOISE = 0           # additional noise, 0 deterministic - 1.0 random noise, ignored for PLMS (i.e. --ddim_eta)
        SCALE = 7.5         # guidance scale, 5 further -> 15 closer to prompt, default 7.5
        SAMPLER = 'PLMS'    # sampler, default 'PLMS' for txt2img which will set NOISE=0, and anything else for DDIM
        STRENGTH = 0.75     # img2img strength, 0 more like original -> 1 less like, default 0.75
        CONFIG = 'configs/stable-diffusion/v1-inference.yaml'
        MODEL = 'models/ldm/stable-diffusion-v1/sd-v1-4.ckpt'

    class IP(SD):           # Inpainting without mask:
        CONFIG = 'models/ldm/inpainting_big/config.yaml'
        MODEL = 'models/ldm/inpainting_big/last.ckpt'

    class CS:               # CLIPSeg `txt2mask` :
        BLACK = 0.2         # threshold <= to which to convert to black, to feather the edges
        WHITE = 0.5         # threshold >= to which to convert to white, to feather the edges
        DILATE = 20         # mask dilate amount in pixels, default 20
        MODEL = 'clipseg/weights/rd64-uni-refined.pth'
        FILENAME = '%y%m%d_%H%M%Sm.png'

    class RE:               # Real-ESRGAN up-scaling, change RealESRGAN_x4plus to RealESRGAN_x4plus_anime_6B for anime
        MODEL = 'realesrgan/weights/RealESRGAN_x4plus.pth' # or 'realesrgan/weights/RealESRGAN_x4plus_anime_6B.pth' only
        FACE_MODEL = 'realesrgan/weights/GFPGANv1.4.pth'
        FILENAME = '%y%m%d_%H%M%Su.png'

    class METADATA:
        SAVE_LEVEL = 1      # save config to PNG metadata, default 0=no metadata, 1=prompt + seed, 2=prompt + seed + SD parameters
        LOAD_FILE = None    # filename of PNG to retrieve config from, default None

# OVERRIDES: EXAMPLES / TEST CASES

    # 1. Stable Diffusion txt2img
    PROMPTS = ['the quick brown fox jumps over the lazy dog']

    # 2. Stable Diffusion txt2img with tiling
    #PROMPTS = ['multicolored parquet flooring']
    #TILING = True

    # 3. Stable Diffusion img2img
    #PROMPTS = ['a dark dreary cyberpunk city with many neon signs']
    #IMAGE = 'data/inpainting_examples/bertrand-gabioud-CpuFzIsHYJ0.png'

    # 4. Stable Diffusion inpainting
    #IMAGE = 'data/inpainting_examples/bench2.png'
    #MASK = 'data/inpainting_examples/bench2_mask.png'

    # 5. Stable Diffusion inpainting with CLIPSeg text2mask
    #IMAGE = 'data/inpainting_examples/overture-creations-5sI6fQgYIuo.png'
    #MASK = ['dog']

    # 6. Stable Diffuion txt+img2img with CLIPSeg txt2mask
    #IMAGE = 'data/inpainting_examples/overture-creations-5sI6fQgYIuo.png'
    #MASK = ['dog']
    #PROMPTS = ['a small tiger cub sitting on a bench']

    # 7. Real-ESRGAN upscaling and GFPGAN face enhancement only
    #IMAGE = 'assets/rick.jpeg'
    #UPSCALE = 2
    #FACE_ENHANCE = True

    # 8. Stable Diffusion with upscaling
    #PROMPTS = ['a black_knight and a red_knight clashing swords in a lush forrest, realistic, 4k, hd']
    #SEED = 588513860
    #UPSCALE = 2

    # 9. Everything
    #PROMPTS = ['a medevial castle with a large tree in front']
    #IMAGE = 'data/inpainting_examples/bertrand-gabioud-CpuFzIsHYJ0.png'
    #SD.STRENGTH = 0.8
    #MASK = ['buildings']
    #UNMASK = ['tree', 'trunk']
    #CS.DILATE = 50
    #UPSCALE = 2
    #FACE_ENHANCE = False

#######################################
# StableDiffusion txt2img
#######################################
from einops import rearrange
from contextlib import nullcontext
from ldm.util import instantiate_from_config
from ldm.models.diffusion.plms import PLMSSampler
from ldm.models.diffusion.ddim import DDIMSampler
from transformers import logging as transformers_logging

class StableDiffusion_Text():
    """Stable Diffusion text to image (M1 MPS)"""

    def __init__(self, device):
        """Initialize Stable Diffusion with the given config and checkpoint files"""
        logging.debug(f'StableDiffusion_Text Init Config={CONFIG.SD.CONFIG}')
        tic = time.perf_counter()
        config = OmegaConf.load(CONFIG.SD.CONFIG)
        model = instantiate_from_config(config.model)
        pl_sd = torch.load(CONFIG.SD.MODEL, map_location='cpu')
        model.load_state_dict(pl_sd['state_dict'], strict=False)
        model.to(device.type)
        model.eval()
        self.model = model
        self.device = device
        logging.debug(f'StableDiffusion_Text Init Time={(time.perf_counter()-tic):.2f}s')

    def setup_sampler(self):
        """Select the sampler, PLMS or DDIM, and some setup if tiling enabled"""
        if CONFIG.SD.SAMPLER == 'PLMS':
            self.sampler = PLMSSampler(self.model)
            CONFIG.SD.NOISE = 0
        else:
            self.sampler = DDIMSampler(self.model)
            CONFIG.SD.SAMPLER = 'DDIM'
        tiling = 'circular' if CONFIG.TILING else 'zeros'
        for m in self.model.modules():
            if isinstance(m, (torch.nn.Conv2d, torch.nn.ConvTranspose2d)):
                if tiling != m.padding_mode: m.padding_mode = tiling
        CONFIG.SEED = seed_everything(CONFIG.SEED)

    def _generate_sample(self, conditioning, unconditional):
        h = CONFIG.HEIGHT // CONFIG.SD.FACTOR
        w = CONFIG.WIDTH // CONFIG.SD.FACTOR
        shape = [4, h, w]
        start = torch.randn([1, 4, h, w], device='cpu').to(torch.device(self.device.type))
        sample, _ = self.sampler.sample(S=CONFIG.SD.STEPS, conditioning=conditioning, batch_size=1, shape=shape, verbose=False, unconditional_guidance_scale=CONFIG.SD.SCALE, unconditional_conditioning=unconditional, eta=CONFIG.SD.NOISE, x_T=start)
        return sample

    def __save_image(self, image):
        image = 255. * rearrange(image.cpu().numpy(), 'c h w -> h w c')
        img = Image.fromarray(image.astype(np.uint8))
        image_file = os.path.join(CONFIG.OUTPUT_DIR, time.strftime(CONFIG.FILENAME))
        img.save(image_file)
        image_size = os.path.getsize(image_file)
        logging.info(f'  ==> Output Image={image_file}, Size={image_size}')
        if image_size < 1000: logging.warning(f'  ❌ Output image file size is too small ❌')
        return image_file

    def create_image(self, positive:str, negative:str)-> str:
        """Create a single image from the prompts"""
        logging.info(f'  Positive=\'{positive}\'')
        if negative: logging.info(f'  Negative=\'{negative}\'')
        with torch.no_grad():
            precision_scope = nullcontext if self.device.type == 'mps' else torch.autocast
            with precision_scope(self.device.type):
                with self.model.ema_scope():
                    uncond = self.model.get_learned_conditioning(negative) if CONFIG.SD.SCALE != 1.0 else None
                    cond = self.model.get_learned_conditioning(positive)
                    sample = self._generate_sample(cond, uncond)
                    sample = self.model.decode_first_stage(sample)
                    images = torch.clamp((sample + 1.0) / 2.0, min=0.0, max=1.0)
        image_file = self.__save_image(images[0])
        METADATA.save(image_file, positive, negative)
        return image_file

    def _create_images(self, positives, negatives, iterations, reset_seed):
        image_files = []
        if negatives:
            negatives.extend([''] * (len(positives)-len(negatives)))
        else:
            negatives = [''] * len(positives)
        for positive, negative in zip(positives, negatives):
            for iteration in range(iterations):
                if reset_seed: CONFIG.SEED = seed_everything(CONFIG.SEED)
                else: CONFIG.SEED += 1
                image_files.append(self.create_image(positive, negative))
        return image_files

    def create_images(self, positives:[str], negatives:[str]=None, iterations=CONFIG.SD.ITERATIONS, reset_seed=True)-> [str]:
        """Create images, repeating a number of times using the same seed - raison detre"""
        logging.info(f'StableDiffusion_Text Sampler={CONFIG.SD.SAMPLER if CONFIG.SD.SAMPLER else "DDIM"}, Size={CONFIG.WIDTH}x{CONFIG.HEIGHT}, Factor={CONFIG.SD.FACTOR}')
        logging.info(f'  Seed={CONFIG.SEED}, Steps={CONFIG.SD.STEPS}, Noise={CONFIG.SD.NOISE}, Scale={CONFIG.SD.SCALE}{", Tiling" if CONFIG.TILING else ""}')
        tic = time.perf_counter()
        image_files = self._create_images(positives, negatives, iterations, reset_seed)
        logging.debug(f'StableDiffusion_Text Run Time={(time.perf_counter()-tic):.2f}s')
        return image_files

#######################################
# StableDiffusion img2img / text+mask2img / img+mask2img
#######################################
from tqdm import tqdm

class DDIMSampler_Mask(DDIMSampler):
    """Override original DDIMSampler to add z_mask and x0 parameters"""

    @torch.no_grad()
    def decode(self, x_latent, cond, t_start, unconditional_guidance_scale=1.0, unconditional_conditioning=None, use_original_steps=False, z_mask = None, x0=None):
        timesteps = np.arange(self.ddpm_num_timesteps) if use_original_steps else self.ddim_timesteps
        timesteps = timesteps[:t_start]
        time_range = np.flip(timesteps)
        total_steps = timesteps.shape[0]
        print(f"Running DDIM Sampling with {total_steps} timesteps")
        iterator = tqdm(time_range, desc='Decoding image', total=total_steps)
        x_dec = x_latent
        for i, step in enumerate(iterator):
            index = total_steps - i - 1
            ts = torch.full((x_latent.shape[0],), step, device=x_latent.device, dtype=torch.long)
            if z_mask is not None and i < total_steps - 2:
                img_orig = self.model.q_sample(x0, ts)
                mask_inv = 1. - z_mask
                x_dec = (img_orig * mask_inv) + (z_mask * x_dec)
            x_dec, _ = self.p_sample_ddim(x_dec, cond, ts, index=index, use_original_steps=use_original_steps,
                unconditional_guidance_scale=unconditional_guidance_scale,
                unconditional_conditioning=unconditional_conditioning)
        return x_dec

class StableDiffusion_Image(StableDiffusion_Text):
    """Stable Diffusion text+image/mask to image (M1 MPS)"""

    def load_image(self, image_file, mask_file=None):
        """Load the image and mask in PNG format, 512x512"""
        self.image_file = image_file
        self.mask_file = mask_file
        image = Image.open(image_file).convert('RGB')
        w, h = image.size
        w, h = map(lambda x: x - x % 32, (w, h))
        w, h = max(w, CONFIG.WIDTH), max(h, CONFIG.HEIGHT)
        image = image.resize((w, h), resample=Image.Resampling.LANCZOS)
        image = np.array(image).astype(np.float32) / 255.0
        image = image[None].transpose(0, 3, 1, 2)
        image = torch.from_numpy(image)
        image = 2. * image - 1.
        image = image.to(self.device.type)
        self.image = image
        if mask_file:
            mask = Image.open(mask_file).convert('L')
            w, h = mask.size
            w, h = map(lambda x: x - x % 32, (w, h))
            w, h = max(w, CONFIG.WIDTH), max(h, CONFIG.HEIGHT)
            mask = mask.resize((w // CONFIG.SD.FACTOR, h // CONFIG.SD.FACTOR), resample=Image.Resampling.LANCZOS)
            mask = np.array(mask).astype(np.float32) / 255.0
            mask = np.tile(mask, (4, 1, 1))
            mask = mask[None].transpose(0, 1, 2, 3)
            mask = torch.from_numpy(mask).to(self.device.type)
            self.mask = mask
        else:
            self.mask = None

    def setup_sampler(self):
        """Select the sampler, DDIM only"""
        CONFIG.SD.SAMPLER = 'DDIM'
        self.sampler = DDIMSampler_Mask(self.model)
        CONFIG.SEED = seed_everything(CONFIG.SEED)

    def _generate_sample(self, conditioning, unconditional):
        self.sampler.make_schedule(ddim_num_steps=CONFIG.SD.STEPS, ddim_eta=CONFIG.SD.NOISE, verbose=False)
        t_enc = int(CONFIG.SD.STRENGTH * CONFIG.SD.STEPS)
        x0 = self.model.encode_first_stage(self.image)
        x0 = self.model.get_first_stage_encoding(x0)
        z_enc = self.sampler.stochastic_encode(x0, torch.tensor([t_enc]).to(self.device.type))
        if self.mask_file:
            random = torch.randn(self.mask.shape, device=self.device)
            z_enc = (self.mask * random) + ((1-self.mask) * z_enc)
        sample = self.sampler.decode(z_enc, conditioning, t_enc, unconditional_guidance_scale=CONFIG.SD.SCALE, unconditional_conditioning=unconditional, z_mask=self.mask, x0=x0)
        return sample

    def create_images(self, positives:[str], negatives:[str]=None, iterations=1, reset_seed=True)-> [str]:
        """Create images, repeating a number of times using the same seed - raison detre"""
        logging.info(f'StableDiffusion_Image Sampler={CONFIG.SD.SAMPLER if CONFIG.SD.SAMPLER else "DDIM"}, Size={CONFIG.WIDTH}x{CONFIG.HEIGHT}, Factor={CONFIG.SD.FACTOR}')
        logging.info(f'  Seed={CONFIG.SEED}, Strength={CONFIG.SD.STRENGTH}, Steps={CONFIG.SD.STEPS}, Noise={CONFIG.SD.NOISE}, Scale={CONFIG.SD.SCALE}{", Tiling" if CONFIG.TILING else ""}')
        logging.info(f'  Image={self.image_file}{f", Mask={self.mask_file}" if self.mask_file else ""}')
        tic = time.perf_counter()
        image_files = super()._create_images(positives, negatives, iterations, reset_seed)
        logging.debug(f'StableDiffusion_Image Run Time={(time.perf_counter()-tic):.2f}s')
        return image_files

#######################################
# StableDiffusion inpaint
#######################################

class StableDiffusion_Erase():
    """Stable Diffusion inpainting (CPU)"""

    def __init__(self, device):
        """Initialize Stable Diffusion Erase with the given config and checkpoint files"""
        logging.debug(f'StableDiffusion Init Config={CONFIG.IP.CONFIG}')
        tic = time.perf_counter()
        config = OmegaConf.load(CONFIG.IP.CONFIG)
        model = instantiate_from_config(config.model)
        pl_sd = torch.load(CONFIG.IP.MODEL, map_location='cpu')
        model.load_state_dict(pl_sd['state_dict'], strict=False)
        model.to(device.type)
        self.device = device
        self.model = model
        logging.debug(f'StableDiffusion Run Time={(time.perf_counter()-tic):.2f}s')

    def setup_sampler(self):
        """Select the sampler, DDIM"""
        CONFIG.IP.SAMPLER = 'DDIM'
        self.sampler = DDIMSampler(self.model)
        CONFIG.SEED = seed_everything(CONFIG.SEED)

    def load_image(self, image_file:str, mask_file:str):
        """Load the image and mask PNG files, max 512x512"""
        self.image_file = image_file
        self.mask_file = mask_file
        image = np.array(Image.open(image_file).convert('RGB'))
        image = image.astype(np.float32)/255.0
        image = image[None].transpose(0,3,1,2)
        image = torch.from_numpy(image)
        mask = np.array(Image.open(mask_file).convert('L'))
        mask = mask.astype(np.float32)/255.0
        mask = mask[None, None]
        mask[mask < 0.5] = 0
        mask[mask >= 0.5] = 1
        mask = torch.from_numpy(mask)
        masked_image = (1-mask) * image
        image = image.to(self.device.type)
        image = 2. * image - 1.
        self.image = image
        mask = mask.to(self.device.type)
        mask = 2. * mask - 1.
        self.mask = mask
        masked_image = masked_image.to(self.device.type)
        masked_image = 2. * masked_image - 1.
        self.masked_image = masked_image

    def __save_image(self, image):
        img = Image.fromarray(image.astype(np.uint8))
        image_file = os.path.join(CONFIG.OUTPUT_DIR, time.strftime(CONFIG.FILENAME))
        img.save(image_file)
        logging.info(f'  ==> Output Image={image_file}')
        return image_file

    def create_image(self)-> str:
        """Create image by erasing masked area - raison detre"""
        logging.info(f'StableDiffusion_Erase Sampler={CONFIG.IP.SAMPLER if CONFIG.IP.SAMPLER else "DDIM"}, Steps={CONFIG.IP.STEPS}')
        if self.device.type == 'cpu': logging.warning(f'  StableDiffusion_Erase running on CPU')
        logging.info(f'  Image={self.image_file}{f", Mask={self.mask_file}" if self.mask_file else ""}')
        tic = time.perf_counter()
        with torch.no_grad():
            with self.model.ema_scope():
                cond = self.model.cond_stage_model.encode(self.masked_image)
                uncond = torch.nn.functional.interpolate(self.mask, size=cond.shape[-2:])
                cond = torch.cat((cond, uncond), dim=1)
                shape = (cond.shape[1]-1,)+cond.shape[2:]
                samples, _ = self.sampler.sample(S=CONFIG.IP.STEPS, conditioning=cond, batch_size=cond.shape[0], shape=shape, verbose=False)
                samples = self.model.decode_first_stage(samples)
                image = torch.clamp((self.image+1.0)/2.0, min=0.0, max=1.0)
                mask = torch.clamp((self.mask+1.0)/2.0, min=0.0, max=1.0)
                predicted_image = torch.clamp((samples+1.0)/2.0, min=0.0, max=1.0)
                inpainted = (1-mask)*image+mask*predicted_image
                inpainted = inpainted.cpu().numpy().transpose(0,2,3,1)[0]*255
        image_file = self.__save_image(inpainted)
        logging.debug(f'StableDiffusion_Erase Run Time={(time.perf_counter()-tic):.2f}s')
        return image_file

#######################################
# CLIPSeg
#######################################
import cv2
from torchvision import transforms
from torchvision.utils import save_image as tv_save_image
from clipseg.models.clipseg import CLIPDensePredT

class CLIPSeg_Mask():
    """CLIPSeg for text to mask (CPU)"""

    def __init__(self, device):
        """Initialize CLIPSeg with the given weights file"""
        logging.debug(f'CLIPSeg_Mask Init Weights={CONFIG.CS.MODEL}')
        tic = time.perf_counter()
        model = CLIPDensePredT(version='ViT-B/16', reduce_dim=64, complex_trans_conv=True)
        pl_sd = torch.load(CONFIG.CS.MODEL, map_location='cpu')
        model.load_state_dict(pl_sd, strict=False)
        model.eval()
        self.device = device
        self.model = model
        logging.debug(f'CLIPSeg_Mask Init Time={(time.perf_counter()-tic):.2f}s')

    def __load_image(self, image_file):
        transform = transforms.Compose([
            transforms.ToTensor(),
            transforms.Normalize(mean=[0.485, 0.456, 0.406], std=[0.229, 0.224, 0.225]),
            transforms.Resize((CONFIG.WIDTH, CONFIG.HEIGHT)) ], )
        image = Image.open(image_file)
        return transform(image).unsqueeze(0)

    def __generate_samples(self, image, prompts):
        length = len(prompts)
        with torch.no_grad():
            preds = self.model(image.repeat(length,1,1,1), prompts)[0]
            mask = sum(torch.sigmoid(preds[i][0]) for i in range(length))
            mask = (mask - mask.min()) / mask.max()
            mask = torch.where(mask >= CONFIG.CS.WHITE, 1., mask)
            mask = torch.where(mask <= CONFIG.CS.BLACK, 0., mask)
        return mask

    def __dilate_mask(self, mask_file, dilate):
        mask = cv2.imread(mask_file, 0)
        kernel = np.ones((dilate, dilate), np.uint8)
        mask = cv2.dilate(mask, kernel, iterations=1)
        cv2.imwrite(mask_file, mask)
        return mask_file

    def create_mask(self, image_file:str, positives:[str], negatives:[str]=None)-> str:
        """Create mask from the prompts, in PNG format - raison detre"""
        logging.info(f'CLIPSeg_Mask Image={image_file}, Size={CONFIG.WIDTH}x{CONFIG.HEIGHT}, Dilate={CONFIG.CS.DILATE}, Black={CONFIG.CS.BLACK}, White={CONFIG.CS.WHITE}')
        logging.info(f'  Positive={positives}')
        tic = time.perf_counter()
        image = self.__load_image(image_file)
        mask = self.__generate_samples(image, positives)
        if negatives:
            logging.info(f'  Negative={negatives}')
            mask -= self.__generate_samples(image, negatives)
        mask_file = os.path.join(CONFIG.OUTPUT_DIR, time.strftime(CONFIG.CS.FILENAME))
        tv_save_image(mask, mask_file)
        if CONFIG.CS.DILATE: self.__dilate_mask(mask_file, CONFIG.CS.DILATE)
        logging.info(f'  ==> Output Mask={mask_file}')
        logging.debug(f'CLIPSeg_Mask Run Time={(time.perf_counter()-tic):.2f}s')
        return mask_file

#######################################
# Real-ESRGAN and GFPGAN upscaling
#######################################
from realesrgan import RealESRGANer
from basicsr.archs.rrdbnet_arch import RRDBNet
from gfpgan import GFPGANer

class RealESRGAN_Upscaler():
    """Real-ESRGAN Upscaler (M1 MPS)"""

    def __init__(self, device):
        """Initialize Real-ESRGAN with the given weights"""
        logging.debug(f'RealESRGAN_Upscaler Init Weights={CONFIG.RE.MODEL}')
        tic = time.perf_counter()
        self.device = device
        if 'RealESRGAN_x4plus_anime_6B' in CONFIG.RE.MODEL:
            self.model = RRDBNet(num_in_ch=3, num_out_ch=3, num_feat=64, num_block=6, num_grow_ch=32, scale=4)
        else:
            self.model = RRDBNet(num_in_ch=3, num_out_ch=3, num_feat=64, num_block=23, num_grow_ch=32, scale=4)
        self.model.to(device.type)
        self.upsampler = RealESRGANer(scale=4, model_path=CONFIG.RE.MODEL, model=self.model, device=device.type)
        logging.debug(f'RealESRGAN Init Time={(time.perf_counter()-tic):.2f}s')

    def load_image(self, image_file:str):
        """Load image, applying MPS fix if required"""
        self.image_file = image_file
        image = cv2.imread(image_file)
        if self.device.type == 'mps':
            image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
            image = image.flatten().reshape((image.shape[2], image.shape[0], image.shape[1])).transpose((1,2,0))
            image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
        self.image = image

    def _save_image(self, image):
        image_file = os.path.join(CONFIG.OUTPUT_DIR, time.strftime(CONFIG.RE.FILENAME))
        cv2.imwrite(image_file, image)
        if CONFIG.PROMPTS and CONFIG.METADATA.SAVE_LEVEL:
            prompts = METADATA.load_str(self.image_file, 'prompt', 'negative', 'Software')
            if 'myByways' in prompts[2]: METADATA.save(image_file, prompts[0], prompts[1])
        logging.info(f'  ==> Output Image={image_file}')
        return image_file

    def upscale_image(self, scale=CONFIG.UPSCALE)->str:
        """Upscale the image - raison detre """
        logging.info(f'RealESRGAN_Upscaler Image={self.image_file}, Scale={scale}')
        tic = time.perf_counter()
        image, _ = self.upsampler.enhance(self.image, outscale=scale)
        image_file = self._save_image(image)
        logging.debug(f'RealESRGAN_Upscaler Run Time={(time.perf_counter()-tic):.2f}s')
        return image_file

class RealESRGAN_FaceUpscaler(RealESRGAN_Upscaler):
    """Real-ESRGAN with GFPGAN Face Enhancement (CPU)"""

    def upscale_image(self, scale=CONFIG.UPSCALE):
        """Upscale the image and enhance faces - raison detre) """
        logging.info(f'RealESRGAN_FaceUpscaler Image={self.image_file}, Scale={scale}, Face Enahnce=True')
        if self.device.type == 'cpu': logging.warning(f'  RealESRGAN_FaceUpscaler running on CPU')
        tic = time.perf_counter()
        face_enhancer = GFPGANer(model_path=CONFIG.RE.FACE_MODEL, upscale=scale, bg_upsampler=self.upsampler, device=self.device)
        _, _, image = face_enhancer.enhance(self.image, has_aligned=False, only_center_face=False, paste_back=True)
        image_file = self._save_image(image)
        logging.debug(f'RealESRGAN_FaceUpscaler Run Time={(time.perf_counter()-tic):.2f}s')
        return image_file

#######################################
# PNG metadata
#######################################
from PIL.PngImagePlugin import PngInfo

class METADATA(CONFIG):
    """Stores Stable Diffusion parameters in PNG iTXt chunk (UTF8)"""

    def __save_list(values):
        if not isinstance(values, list): return values
        masks = [mask.replace('|', ' ') for mask in CONFIG.MASK]
        return '|'+'|'.join(masks)+'|'

    def __load_list(values):
        if values[0] != '|': return values
        values = values.split('|')
        return [value for value in values if value]

    def save(image_file, prompt, negative=None, level=CONFIG.METADATA.SAVE_LEVEL):
        metadata = PngInfo()
        metadata.add_itxt('Software', 'Simple-SD v1.0 copyright (c) 2022 C.Y. Wong, myByways.com')
        metadata.add_itxt('Author', 'Stable Diffusion v1')
        if level >= 1:
            metadata.add_itxt('prompt', prompt)
            if negative: metadata.add_itxt('negative', negative)
            metadata.add_itxt('seed', str(CONFIG.SEED))
            if CONFIG.TILING: metadata.add_itxt('tiling', str(1))
        if level >= 2:
            metadata.add_itxt('steps', str(CONFIG.SD.STEPS))
            metadata.add_itxt('noise', str(CONFIG.SD.NOISE))
            metadata.add_itxt('scale',str(CONFIG.SD.SCALE))
            if CONFIG.SD.SAMPLER: metadata.add_itxt('sampler', CONFIG.SD.SAMPLER)
            if CONFIG.UPSCALE: metadata.add_itxt('upscale', str(CONFIG.UPSCALE))
        if level >= 3:
            if CONFIG.IMAGE:
                metadata.add_itxt('image', CONFIG.IMAGE)
                metadata.add_itxt('strength', str(CONFIG.SD.STRENGTH))
            if CONFIG.MASK: metadata.add_itxt('mask', METADATA.__save_list(CONFIG.MASK))
            if CONFIG.UNMASK: metadata.add_itxt('unmask', METADATA.__save_list(CONFIG.UNMASK))
        if image_file.lower().endswith('.png'):
            image = Image.open(image_file)
            image.save(image_file, pnginfo=metadata)
        return metadata

    def load_str(image_file, *keys):
        values = []
        image = Image.open(image_file)
        for key in keys:
            values.append(image.text.get(key, None))
        return values

    def load(image_file=CONFIG.METADATA.LOAD_FILE):
        logging.debug(f'METADATA Loading from PNG File={image_file}')
        if not image_file.lower().endswith('.png'): return None
        image = Image.open(image_file)
        CONFIG.TILING = False
        CONFIG.SD.SAMPLER = ''
        for k, v in image.text.items():
            if k == 'prompt': CONFIG.PROMPTS = [v]
            elif k == 'negative': CONFIG.NEGATIVES = [v]
            elif k == 'seed': CONFIG.SEED = int(v) if v != 'None' else None
            elif k == 'tiling': CONFIG.TILING = True
            elif k == 'steps': CONFIG.SD.STEPS = int(v)
            elif k == 'noise': CONFIG.SD.NOISE = float(v)
            elif k == 'scale': CONFIG.SD.SCALE = float(v)
            elif k == 'sampler': CONFIG.SD.SAMPLER = v
            elif k == 'upscale': CONFIG.UPSCALE = float(v)
            elif k == 'image': CONFIG.IMAGE = v
            elif k == 'strength': CONFIG.SD.STRENGTH = float(v)
            elif k == 'mask': CONFIG.MASK = METADATA.__load_list(v)
            elif k == 'unmask': CONFIG.UNMASK = METADATA.__load_list(v)
        return image.text

    def __str__(self):
        return f'Prompts = {CONFIG.PROMPTS}\nNegatives = {CONFIG.NEGATIVES}\nTiling = {CONFIG.TILING}\n' \
            f'Image = {CONFIG.IMAGE}\nMask = {CONFIG.MASK}\nUnmask = {CONFIG.UNMASK}\n' \
            f'Upscale = {CONFIG.UPSCALE}\nFace_Enhance = {CONFIG.FACE_ENHANCE}\nSeed = {CONFIG.SEED}\n' \
            f'Width = {CONFIG.WIDTH}\nHeight = {CONFIG.HEIGHT}\nFilename = {CONFIG.FILENAME}\n' \
            f'SD.Factor = {CONFIG.SD.FACTOR}\nSD.Steps = {CONFIG.SD.STEPS}\nSD.Noise = {CONFIG.SD.NOISE}\n' \
            f'SD.Scale = {CONFIG.SD.SCALE}\nSD.Strength = {CONFIG.SD.STRENGTH}\nSD.Config = {CONFIG.SD.CONFIG}\n' \
            f'SD.Iterations = {CONFIG.SD.ITERATIONS}\nIP.Config = {CONFIG.IP.CONFIG}\nIP.Model = {CONFIG.IP.MODEL}\n' \
            f'CS.Black = {CONFIG.CS.BLACK}\nCS.White = {CONFIG.CS.WHITE}\nCS.Dilage = {CONFIG.CS.DILATE}\n' \
            f'CS.Model = {CONFIG.CS.MODEL}\nCS.Filename = {CONFIG.CS.FILENAME}\nRE.Model = {CONFIG.RE.MODEL}\n' \
            f'RE.Face_Model = {CONFIG.RE.FACE_MODEL}\nRE.Filename = {CONFIG.RE.FILENAME}'

def setup_logging():
    class ConsoleLog_Formatter(logging.Formatter):
        formats = {
            logging.DEBUG:    logging.Formatter('🔵 \x1b[34m%(message)s\x1b[0m'),
            logging.INFO:     logging.Formatter('🟢 \x1b[32m%(message)s\x1b[0m'),
            logging.WARNING:  logging.Formatter('🟡 \x1b[33;1m%(message)s\x1b[0m'),
            logging.ERROR:    logging.Formatter('❌ \x1b[31;1m%(message)s\x1b[0m'),
            logging.CRITICAL: logging.Formatter('❌ \x1b[31;1m%(message)s\x1b[0m')
        }
        def format(self, log_record):
            return self.formats[log_record.levelno].format(log_record)

    class FileLog_InfoFilter(object):
        def filter(self, log_record):
            return log_record.levelno <= logging.INFO

    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    if len(logger.root.handlers) == 0:
        logger.addHandler(logging.StreamHandler())
    logger.handlers[0].setLevel(logging.DEBUG)
    logger.handlers[0].setFormatter(ConsoleLog_Formatter())
    file_log = logging.FileHandler(os.path.join(CONFIG.OUTPUT_DIR, CONFIG.HISTORY_FILE))
    file_log.setLevel(logging.INFO)
    file_log.addFilter(FileLog_InfoFilter())
    file_log.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(file_log)
    warnings.filterwarnings('ignore')
    logging.getLogger('PIL.PngImagePlugin').setLevel(logging.ERROR)
    logging.getLogger('pytorch_lightning').setLevel(logging.ERROR)
    transformers_logging.set_verbosity_error()

#######################################
# Main
#######################################

def main():
    print('\x1b[44m🙈🙉🙊 \x1b[1;34m\x1b[4mmyByways Simple Stable Diffusion v1.0 (1 Oct 2022)\x1b[0m\x1b[44m 🙈🙉🙊\x1b[0m')
    setup_logging()
    os.makedirs(CONFIG.OUTPUT_DIR, exist_ok=True)

    if CONFIG.METADATA.LOAD_FILE:
        METADATA.load(CONFIG.METADATA.LOAD_FILE)

    if torch.backends.mps.is_available():
        device = torch.device('mps')
    elif torch.cuda.is_available():
        device = torch.device('cuda')
    else:
        device = torch.device('cpu')

    mask_file = CONFIG.MASK
    if CONFIG.MASK and isinstance(CONFIG.MASK, list):
        cs = CLIPSeg_Mask(device)
        mask_file = cs.create_mask(CONFIG.IMAGE, CONFIG.MASK, CONFIG.UNMASK)

    image_files = None
    if CONFIG.PROMPTS:
        if CONFIG.IMAGE:
            sd = StableDiffusion_Image(device)
            sd.load_image(CONFIG.IMAGE, mask_file)
        else:
            sd = StableDiffusion_Text(device)
        sd.setup_sampler()
        image_files = sd.create_images(CONFIG.PROMPTS, CONFIG.NEGATIVES)
    elif CONFIG.IMAGE and mask_file:
        sd = StableDiffusion_Erase(torch.device('cpu'))
        sd.setup_sampler()
        sd.load_image(CONFIG.IMAGE, mask_file)
        image_files = [sd.create_image()]
    else:
        image_files = [CONFIG.IMAGE]

    upscaled_files = []
    if CONFIG.UPSCALE and image_files:
        if CONFIG.FACE_ENHANCE:
            re = RealESRGAN_FaceUpscaler(torch.device('cpu'))
        else:
            re = RealESRGAN_Upscaler(device)
        for image_file in image_files:
            re.load_image(image_file)
            upscaled_files.append(re.upscale_image(CONFIG.UPSCALE))

    if mask_file != CONFIG.MASK: logging.debug(f'😎 Mask     = [{mask_file}]')
    if image_files:              logging.debug(f'😎 Images   = {image_files}')
    if upscaled_files:           logging.debug(f'😎 Upscaled = {upscaled_files}')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        logging.error('User abort, goodbye')
    except FileNotFoundError as e:
        logging.error(f'File not found {e}')
