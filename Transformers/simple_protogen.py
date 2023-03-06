from diffusers import StableDiffusionPipeline, DPMSolverMultistepScheduler
import torch

prompt = (
"modelshoot style, (extremely detailed CG unity 8k wallpaper), full shot body photo of the most beautiful artwork in the world, "
"english medieval witch, black silk vale, pale skin, black silk robe, black cat, necromancy magic, medieval era, "
"photorealistic painting by Ed Blinkey, Atey Ghailan, Studio Ghibli, by Jeremy Mann, Greg Manchess, Antonio Moro, trending on ArtStation, "
"trending on CGSociety, Intricate, High Detail, Sharp focus, dramatic, photorealistic painting art by midjourney and greg rutkowski"
)

model_id = "darkstorm2150/Protogen_x3.4_Official_Release"
pipe = StableDiffusionPipeline.from_pretrained(model_id, torch_dtype=torch.float16)
pipe.scheduler = DPMSolverMultistepScheduler.from_config(pipe.scheduler.config)
pipe = pipe.to("cuda")

image = pipe(prompt, num_inference_steps=25).images[0]

image.save("./result.jpg")
