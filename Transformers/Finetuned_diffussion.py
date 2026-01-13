#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  finetuned_diffusion.py
#
#  Copyright 2023 Diego Martinez Gutierrez <diego.martinez@ehu.eus>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
#
# ---------------------------
# Importación de los módulos
# ---------------------------
from diffusers import AutoencoderKL, UNet2DConditionModel, StableDiffusionPipeline, StableDiffusionImg2ImgPipeline, DPMSolverMultistepScheduler
import gradio as gr
import torch
from PIL import Image
import utils
import datetime
import time
import psutil
import random


start_time = time.time()
is_colab = utils.is_google_colab()
state = None
current_steps = 25

class Model:
    def __init__(self, name, path="", prefix=""):
        self.name = name
        self.path = path
        self.prefix = prefix
        self.pipe_t2i = None
        self.pipe_i2i = None

models = [
     Model("Arcane", "nitrosocke/Arcane-Diffusion", "arcane style "),
     Model("Dreamlike Diffusion 1.0", "dreamlike-art/dreamlike-diffusion-1.0", "dreamlikeart "),
     Model("Archer", "nitrosocke/archer-diffusion", "archer style "),
     Model("Anything V3 Better-Vae", "Linaqruf/anything-v3-better-vae", ""),
     Model("Anything V4", "andite/anything-v4.0", ""),
     Model("Anything midjourney", "Joeythemonster/anything-midjourney-v-4-1", ""),
     Model("Anything V4.5", "ckpt/anything-v4.5", ""),
     Model( "Pastel Mix","andite/pastel-mix", ""),
     Model("Yohan","andite/yohan-diffusion", ""),
     Model("cutesexyrobutts","andite/cutesexyrobutts-diffusion",""),
     Model("Modern Disney", "nitrosocke/mo-di-diffusion", "modern disney style "),
     Model("Classic Disney", "nitrosocke/classic-anim-diffusion", "classic disney style "),
     Model("Loving Vincent (Van Gogh)", "dallinmackay/Van-Gogh-diffusion", "lvngvncnt "),
     Model("Wavyfusion", "wavymulder/wavyfusion", "wa-vy style "),
     Model("Analog Diffusion", "wavymulder/Analog-Diffusion", "analog style "),
     Model("Redshift renderer (Cinema4D)", "nitrosocke/redshift-diffusion", "redshift style "),
     Model("Midjourney v4 style", "prompthero/midjourney-v4-diffusion", "mdjrny-v4 style "),
     Model("Waifu", "hakurei/waifu-diffusion"),
     Model("Cyberpunk Anime", "DGSpitzer/Cyberpunk-Anime-Diffusion", "dgs illustration style "),
     Model("Elden Ring", "nitrosocke/elden-ring-diffusion", "elden ring style "),
     Model("TrinArt v2", "naclbit/trinart_stable_diffusion_v2"),
     Model("Spider-Verse", "nitrosocke/spider-verse-diffusion", "spiderverse style "),
     Model("Balloon Art", "Fictiverse/Stable_Diffusion_BalloonArt_Model", "BalloonArt "),
     Model("Tron Legacy", "dallinmackay/Tron-Legacy-diffusion", "trnlgcy "),
     Model("Pokémon", "lambdalabs/sd-pokemon-diffusers"),
     Model("Pony Diffusion", "AstraliteHeart/pony-diffusion"),
     Model("Robo Diffusion", "nousr/robo-diffusion"),
     Model("Epic Diffusion", "johnslegers/epic-diffusion")
  ]

custom_model = None
if is_colab:
  models.insert(0, Model("Custom model"))
  custom_model = models[0]

last_mode = "txt2img"
current_model = models[1] if is_colab else models[0]
current_model_path = current_model.path

if is_colab:
  pipe = StableDiffusionPipeline.from_pretrained(
      current_model.path,
      torch_dtype=torch.float16,
      scheduler=DPMSolverMultistepScheduler.from_pretrained(current_model.path, subfolder="scheduler"),
      safety_checker=lambda images, clip_input: (images, False)
      )

else:
  pipe = StableDiffusionPipeline.from_pretrained(
      current_model.path,
      torch_dtype=torch.float16,
      scheduler=DPMSolverMultistepScheduler.from_pretrained(current_model.path, subfolder="scheduler")
      )

if torch.cuda.is_available():
  pipe = pipe.to("cuda")
  pipe.enable_xformers_memory_efficient_attention()

device = "GPU 🔥" if torch.cuda.is_available() else "CPU 🥶"

def error_str(error, title="Error"):
    return f"""#### {title}
            {error}"""  if error else ""

def update_state(new_state):
  global state
  state = new_state

def update_state_info(old_state):
  if state and state != old_state:
    return gr.update(value=state)

def custom_model_changed(path):
  models[0].path = path
  global current_model
  current_model = models[0]

def on_model_change(model_name):

  prefix = "Enter prompt. \"" + next((m.prefix for m in models if m.name == model_name), None) + "\" is prefixed automatically" if model_name != models[0].name else "Don't forget to use the custom model prefix in the prompt!"

  return gr.update(visible = model_name == models[0].name), gr.update(placeholder=prefix)

def on_steps_change(steps):
  global current_steps
  current_steps = steps

def pipe_callback(step: int, timestep: int, latents: torch.FloatTensor):
    update_state(f"{step}/{current_steps} steps")#\nTime left, sec: {timestep/100:.0f}")

def inference(model_name, prompt, guidance, steps, n_images=1, width=512, height=512, seed=0, img=None, strength=0.5, neg_prompt=""):

  update_state(" ")

  print(psutil.virtual_memory()) # print memory usage

  global current_model
  for model in models:
    if model.name == model_name:
      current_model = model
      model_path = current_model.path

  # generator = torch.Generator('cuda').manual_seed(seed) if seed != 0 else None
  if seed == 0:
    seed = random.randint(0, 2147483647)

  generator = torch.Generator('cuda').manual_seed(seed)

  try:
    if img is not None:
      return img_to_img(model_path, prompt, n_images, neg_prompt, img, strength, guidance, steps, width, height, generator, seed), f"Done. Seed: {seed}"
    else:
      return txt_to_img(model_path, prompt, n_images, neg_prompt, guidance, steps, width, height, generator, seed), f"Done. Seed: {seed}"
  except Exception as e:
    return None, error_str(e)

def txt_to_img(model_path, prompt, n_images, neg_prompt, guidance, steps, width, height, generator, seed):

    print(f"{datetime.datetime.now()} txt_to_img, model: {current_model.name}")

    global last_mode
    global pipe
    global current_model_path
    if model_path != current_model_path or last_mode != "txt2img":
        current_model_path = model_path

        update_state(f"Loading {current_model.name} text-to-image model...")

        if is_colab or current_model == custom_model:
          pipe = StableDiffusionPipeline.from_pretrained(
              current_model_path,
              torch_dtype=torch.float16,
              scheduler=DPMSolverMultistepScheduler.from_pretrained(current_model.path, subfolder="scheduler"),
              safety_checker=lambda images, clip_input: (images, False)
              )
        else:
          pipe = StableDiffusionPipeline.from_pretrained(
              current_model_path,
              torch_dtype=torch.float16,
              scheduler=DPMSolverMultistepScheduler.from_pretrained(current_model.path, subfolder="scheduler")
              )
          # pipe = pipe.to("cpu")
          # pipe = current_model.pipe_t2i

        if torch.cuda.is_available():
          pipe = pipe.to("cuda")
          pipe.enable_xformers_memory_efficient_attention()
        last_mode = "txt2img"

    prompt = current_model.prefix + prompt
    result = pipe(
      prompt,
      negative_prompt = neg_prompt,
      num_images_per_prompt=n_images,
      num_inference_steps = int(steps),
      guidance_scale = guidance,
      width = width,
      height = height,
      generator = generator,
      callback=pipe_callback)

    # update_state(f"Done. Seed: {seed}")

    return replace_nsfw_images(result)

def img_to_img(model_path, prompt, n_images, neg_prompt, img, strength, guidance, steps, width, height, generator, seed):

    print(f"{datetime.datetime.now()} img_to_img, model: {model_path}")

    global last_mode
    global pipe
    global current_model_path
    if model_path != current_model_path or last_mode != "img2img":
        current_model_path = model_path

        update_state(f"Loading {current_model.name} image-to-image model...")

        if is_colab or current_model == custom_model:
          pipe = StableDiffusionImg2ImgPipeline.from_pretrained(
              current_model_path,
              torch_dtype=torch.float16,
              scheduler=DPMSolverMultistepScheduler.from_pretrained(current_model.path, subfolder="scheduler"),
              safety_checker=lambda images, clip_input: (images, False)
              )
        else:
          pipe = StableDiffusionImg2ImgPipeline.from_pretrained(
              current_model_path,
              torch_dtype=torch.float16,
              scheduler=DPMSolverMultistepScheduler.from_pretrained(current_model.path, subfolder="scheduler")
              )
          # pipe = pipe.to("cpu")
          # pipe = current_model.pipe_i2i

        if torch.cuda.is_available():
          pipe = pipe.to("cuda")
          pipe.enable_xformers_memory_efficient_attention()
        last_mode = "img2img"

    prompt = current_model.prefix + prompt
    ratio = min(height / img.height, width / img.width)
    img = img.resize((int(img.width * ratio), int(img.height * ratio)), Image.LANCZOS)
    result = pipe(
        prompt,
        negative_prompt = neg_prompt,
        num_images_per_prompt=n_images,
        image = img,
        num_inference_steps = int(steps),
        strength = strength,
        guidance_scale = guidance,
        # width = width,
        # height = height,
        generator = generator,
        callback=pipe_callback)

    # update_state(f"Done. Seed: {seed}")

    return replace_nsfw_images(result)

def replace_nsfw_images(results):

    if is_colab:
      return results.images

    for i in range(len(results.images)):
      if results.nsfw_content_detected[i]:
        results.images[i] = Image.open("nsfw.png")
    return results.images

# css = """.finetuned-diffusion-div div{display:inline-flex;align-items:center;gap:.8rem;font-size:1.75rem}.finetuned-diffusion-div div h1{font-weight:900;margin-bottom:7px}.finetuned-diffusion-div p{margin-bottom:10px;font-size:94%}a{text-decoration:underline}.tabs{margin-top:0;margin-bottom:0}#gallery{min-height:20rem}
# """
with gr.Blocks(css="style.css") as demo:
    gr.HTML(
        f"""
            <div class="finetuned-diffusion-div">
              <div>
                <h1>Finetuned Diffusion</h1>
              </div>
              <p>
               Demo for multiple fine-tuned Stable Diffusion models, trained on different styles: <br>
               <a href="https://huggingface.co/nitrosocke/Arcane-Diffusion">Arcane</a>, <a href="https://huggingface.co/nitrosocke/archer-diffusion">Archer</a>, <a href="https://huggingface.co/nitrosocke/elden-ring-diffusion">Elden Ring</a>, <a href="https://huggingface.co/nitrosocke/spider-verse-diffusion">Spider-Verse</a>, <a href="https://huggingface.co/nitrosocke/mo-di-diffusion">Modern Disney</a>, <a href="https://huggingface.co/nitrosocke/classic-anim-diffusion">Classic Disney</a>, <a href="https://huggingface.co/dallinmackay/Van-Gogh-diffusion">Loving Vincent (Van Gogh)</a>, <a href="https://huggingface.co/nitrosocke/redshift-diffusion">Redshift renderer (Cinema4D)</a>, <a href="https://huggingface.co/prompthero/midjourney-v4-diffusion">Midjourney v4 style</a>, <a href="https://huggingface.co/hakurei/waifu-diffusion">Waifu</a>, <a href="https://huggingface.co/lambdalabs/sd-pokemon-diffusers">Pokémon</a>, <a href="https://huggingface.co/AstraliteHeart/pony-diffusion">Pony Diffusion</a>, <a href="https://huggingface.co/nousr/robo-diffusion">Robo Diffusion</a>, <a href="https://huggingface.co/DGSpitzer/Cyberpunk-Anime-Diffusion">Cyberpunk Anime</a>, <a href="https://huggingface.co/dallinmackay/Tron-Legacy-diffusion">Tron Legacy</a>, <a href="https://huggingface.co/Fictiverse/Stable_Diffusion_BalloonArt_Model">Balloon Art</a> + in colab notebook you can load any other Diffusers 🧨 SD model hosted on HuggingFace 🤗.
              </p>
              <p>You can skip the queue and load custom models in the colab: <a href="https://colab.research.google.com/gist/qunash/42112fb104509c24fd3aa6d1c11dd6e0/copy-of-fine-tuned-diffusion-gradio.ipynb"><img data-canonical-src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab" src="https://camo.githubusercontent.com/84f0493939e0c4de4e6dbe113251b4bfb5353e57134ffd9fcab6b8714514d4d1/68747470733a2f2f636f6c61622e72657365617263682e676f6f676c652e636f6d2f6173736574732f636f6c61622d62616467652e737667"></a></p>
               Running on <b>{device}</b>{(" in a <b>Google Colab</b>." if is_colab else "")}
              </p>
              <p>You can also duplicate this space and upgrade to gpu by going to settings:<br>
              <a style="display:inline-block" href="https://huggingface.co/spaces/anzorq/finetuned_diffusion?duplicate=true"><img src="https://img.shields.io/badge/-Duplicate%20Space-blue?labelColor=white&style=flat&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAAAXNSR0IArs4c6QAAAP5JREFUOE+lk7FqAkEURY+ltunEgFXS2sZGIbXfEPdLlnxJyDdYB62sbbUKpLbVNhyYFzbrrA74YJlh9r079973psed0cvUD4A+4HoCjsA85X0Dfn/RBLBgBDxnQPfAEJgBY+A9gALA4tcbamSzS4xq4FOQAJgCDwV2CPKV8tZAJcAjMMkUe1vX+U+SMhfAJEHasQIWmXNN3abzDwHUrgcRGmYcgKe0bxrblHEB4E/pndMazNpSZGcsZdBlYJcEL9Afo75molJyM2FxmPgmgPqlWNLGfwZGG6UiyEvLzHYDmoPkDDiNm9JR9uboiONcBXrpY1qmgs21x1QwyZcpvxt9NS09PlsPAAAAAElFTkSuQmCC&logoWidth=14" alt="Duplicate Space"></a></p>
            </div>
        """
    )
    with gr.Row():

        with gr.Column(scale=55):
          with gr.Group():
              model_name = gr.Dropdown(label="Model", choices=[m.name for m in models], value=current_model.name)
              with gr.Box(visible=False) as custom_model_group:
                custom_model_path = gr.Textbox(label="Custom model path", placeholder="Path to model, e.g. nitrosocke/Arcane-Diffusion", interactive=True)
                gr.HTML("<div><font size='2'>Custom models have to be downloaded first, so give it some time.</font></div>")

              with gr.Row():
                prompt = gr.Textbox(label="Prompt", show_label=False, max_lines=2,placeholder="Enter prompt. Style applied automatically").style(container=False)
                generate = gr.Button(value="Generate").style(rounded=(False, True, True, False))


              # image_out = gr.Image(height=512)
              gallery = gr.Gallery(label="Generated images", show_label=False, elem_id="gallery").style(grid=[2], height="auto")

          state_info = gr.Textbox(label="State", show_label=False, max_lines=2).style(container=False)
          error_output = gr.Markdown()

        with gr.Column(scale=45):
          with gr.Tab("Options"):
            with gr.Group():
              neg_prompt = gr.Textbox(label="Negative prompt", placeholder="What to exclude from the image")

              n_images = gr.Slider(label="Images", value=1, minimum=1, maximum=4, step=1)

              with gr.Row():
                guidance = gr.Slider(label="Guidance scale", value=7.5, maximum=15)
                steps = gr.Slider(label="Steps", value=current_steps, minimum=2, maximum=75, step=1)

              with gr.Row():
                width = gr.Slider(label="Width", value=512, minimum=64, maximum=1024, step=8)
                height = gr.Slider(label="Height", value=512, minimum=64, maximum=1024, step=8)

              seed = gr.Slider(0, 2147483647, label='Seed (0 = random)', value=0, step=1)

          with gr.Tab("Image to image"):
              with gr.Group():
                image = gr.Image(label="Image", height=256, tool="editor", type="pil")
                strength = gr.Slider(label="Transformation strength", minimum=0, maximum=1, step=0.01, value=0.5)

    if is_colab:
        model_name.change(on_model_change, inputs=model_name, outputs=[custom_model_group, prompt], queue=False)
        custom_model_path.change(custom_model_changed, inputs=custom_model_path, outputs=None)
    # n_images.change(lambda n: gr.Gallery().style(grid=[2 if n > 1 else 1], height="auto"), inputs=n_images, outputs=gallery)
    steps.change(on_steps_change, inputs=[steps], outputs=[], queue=False)

    inputs = [model_name, prompt, guidance, steps, n_images, width, height, seed, image, strength, neg_prompt]
    outputs = [gallery, error_output]
    prompt.submit(inference, inputs=inputs, outputs=outputs)
    generate.click(inference, inputs=inputs, outputs=outputs)

    ex = gr.Examples([
        [models[7].name, "tiny cute and adorable kitten adventurer dressed in a warm overcoat with survival gear on a winters day", 7.5, 25],
        [models[4].name, "portrait of dwayne johnson", 7.0, 35],
        [models[5].name, "portrait of a beautiful alyx vance half life", 10, 25],
        [models[6].name, "Aloy from Horizon: Zero Dawn, half body portrait, smooth, detailed armor, beautiful face, illustration", 7.0, 30],
        [models[5].name, "fantasy portrait painting, digital art", 4.0, 20],
    ], inputs=[model_name, prompt, guidance, steps], outputs=outputs, fn=inference, cache_examples=False)

    gr.HTML("""
    <div style="border-top: 1px solid #303030;">
      <br>
      <p>Models by <a href="https://huggingface.co/nitrosocke">@nitrosocke</a>, <a href="https://twitter.com/haruu1367">@haruu1367</a>, <a href="https://twitter.com/DGSpitzer">@Helixngc7293</a>, <a href="https://twitter.com/dal_mack">@dal_mack</a>, <a href="https://twitter.com/prompthero">@prompthero</a> and others. ❤️</p>
      <p>This space uses the <a href="https://github.com/LuChengTHU/dpm-solver">DPM-Solver++</a> sampler by <a href="https://arxiv.org/abs/2206.00927">Cheng Lu, et al.</a>.</p>
      <p>Space by:<br>
      <a href="https://twitter.com/hahahahohohe"><img src="https://img.shields.io/twitter/follow/hahahahohohe?label=%40anzorq&style=social" alt="Twitter Follow"></a><br>
      <a href="https://github.com/qunash"><img alt="GitHub followers" src="https://img.shields.io/github/followers/qunash?style=social" alt="Github Follow"></a></p><br><br>
      <a href="https://www.buymeacoffee.com/anzorq" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/v2/default-yellow.png" alt="Buy Me A Coffee" style="height: 45px !important;width: 162px !important;" ></a><br><br>
      <p><img src="https://visitor-badge.glitch.me/badge?page_id=anzorq.finetuned_diffusion" alt="visitors"></p>
    </div>
    """)

    demo.load(update_state_info, inputs=state_info, outputs=state_info, every=0.5, show_progress=False)

print(f"Space built in {time.time() - start_time:.2f} seconds")

# if not is_colab:
demo.queue(concurrency_count=1)
demo.launch(debug=is_colab, share=is_colab)
