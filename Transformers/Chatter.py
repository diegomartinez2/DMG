#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  untitled.py
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

from transformers import AutoModelForCausalLM, AutoTokenizer
import torch
import time

# -------
# Clases
# -------
class NombredeClase(object):
    """docstring for NombredeClase."""

    def __init__(self, model_name):
        super(NombredeClase, self).__init__()
        self.model_name = model_name
        # model_name = "microsoft/DialoGPT-large"
        # model_name = "microsoft/DialoGPT-medium"
        # model_name = "microsoft/DialoGPT-small"
        # model_name = "mattallio/Archivist-medium-dialoGPT"
        model_name = "EZSNoVa/DialogGPT-medium-NoVa"
        # model_name = murtawh/DialoGPT-medium-brainarchangel
        # model_name = TheHappyDrone/DialoGPT-medium-Nexus-Nova-turing-v2
        # model_name = allenai/cosmo-xl
        #? model_name = EleutherAI/gpt-neox-20b #don't know it this works
        #? model_name = EleutherAI/gpt-j-6B
        #? model_name = EleutherAI/pythia-6.9b-deduped
        tokenizer = AutoTokenizer.from_pretrained(model_name)
        model = AutoModelForCausalLM.from_pretrained(model_name)

    def run(self, time_limit=300):
        #for step in range(5):
        # while True:
        finished = False
        startTime = time.time()
        while not finished:
            # take user input
            text = input(">> You:")
            # encode the input and add end of string token
            input_ids = tokenizer.encode(text + tokenizer.eos_token, return_tensors="pt")
            # concatenate new user input with chat history (if there is)
            bot_input_ids = torch.cat([chat_history_ids, input_ids], dim=-1) if step > 0 else input_ids
            # generate a bot response
            chat_history_ids = model.generate(
                bot_input_ids,
                max_length=1000,
                do_sample=True,
                top_p=0.95,
                top_k=0,
                temperature=0.75,
                pad_token_id=tokenizer.eos_token_id
            )
            #print the output
            output = tokenizer.decode(chat_history_ids[:, bot_input_ids.shape[-1]:][0], skip_special_tokens=True)
            print(f"DialoGPT: {output}")
            # if condition():
            #     break
            executionTime = (time.time() - startTime)
            # if (executionTime-time_limit > 0):
                # evaluate_end_condition = True
            # elif:
            #     evaluate_end_condition = False
            # finished = evaluate_end_condition
            if (executionTime-time_limit > 0):
                finished = True
        return 0

# ----------
# Funciones
# ----------
def Testing_chat(arg):
    # model_name = "microsoft/DialoGPT-large"
    model_name = "microsoft/DialoGPT-medium"
    # model_name = "microsoft/DialoGPT-small"
    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModelForCausalLM.from_pretrained(model_name)

    # chatting 5 times with greedy search
    for step in range(arg):
        # take user input
        text = input(">> You:")
        # encode the input and add end of string token
        input_ids = tokenizer.encode(text + tokenizer.eos_token, return_tensors="pt")
        # concatenate new user input with chat history (if there is)
        bot_input_ids = torch.cat([chat_history_ids, input_ids], dim=-1) if step > 0 else input_ids
        # generate a bot response
        chat_history_ids = model.generate(
            bot_input_ids,
            max_length=1000,
            pad_token_id=tokenizer.eos_token_id,
        )
        #print the output
        output = tokenizer.decode(chat_history_ids[:, bot_input_ids.shape[-1]:][0], skip_special_tokens=True)
        print(f"DialoGPT: {output}")

    # chatting 5 times with beam search
    for step in range(arg):
        # take user input
        text = input(">> You:")
        # encode the input and add end of string token
        input_ids = tokenizer.encode(text + tokenizer.eos_token, return_tensors="pt")
        # concatenate new user input with chat history (if there is)
        bot_input_ids = torch.cat([chat_history_ids, input_ids], dim=-1) if step > 0 else input_ids
        # generate a bot response
        chat_history_ids = model.generate(
            bot_input_ids,
            max_length=1000,
            num_beams=3,
            early_stopping=True,
            pad_token_id=tokenizer.eos_token_id
        )
        #print the output
        output = tokenizer.decode(chat_history_ids[:, bot_input_ids.shape[-1]:][0], skip_special_tokens=True)
        print(f"DialoGPT: {output}")

    # chatting 5 times with sampling
    for step in range(arg):
        # take user input
        text = input(">> You:")
        # encode the input and add end of string token
        input_ids = tokenizer.encode(text + tokenizer.eos_token, return_tensors="pt")
        # concatenate new user input with chat history (if there is)
        bot_input_ids = torch.cat([chat_history_ids, input_ids], dim=-1) if step > 0 else input_ids
        # generate a bot response
        chat_history_ids = model.generate(
            bot_input_ids,
            max_length=1000,
            do_sample=True,
            top_k=0,
            pad_token_id=tokenizer.eos_token_id
        )
        #print the output
        output = tokenizer.decode(chat_history_ids[:, bot_input_ids.shape[-1]:][0], skip_special_tokens=True)
        print(f"DialoGPT: {output}")

    # chatting 5 times with Top K sampling & tweaking temperature
    for step in range(arg):
        # take user input
        text = input(">> You:")
        # encode the input and add end of string token
        input_ids = tokenizer.encode(text + tokenizer.eos_token, return_tensors="pt")
        # concatenate new user input with chat history (if there is)
        bot_input_ids = torch.cat([chat_history_ids, input_ids], dim=-1) if step > 0 else input_ids
        # generate a bot response
        chat_history_ids = model.generate(
            bot_input_ids,
            max_length=1000,
            do_sample=True,
            top_k=100,
            temperature=0.75,
            pad_token_id=tokenizer.eos_token_id
        )
        #print the output
        output = tokenizer.decode(chat_history_ids[:, bot_input_ids.shape[-1]:][0], skip_special_tokens=True)
        print(f"DialoGPT: {output}")

    # chatting 5 times with nucleus sampling & tweaking temperature
    for step in range(arg):
        # take user input
        text = input(">> You:")
        # encode the input and add end of string token
        input_ids = tokenizer.encode(text + tokenizer.eos_token, return_tensors="pt")
        # concatenate new user input with chat history (if there is)
        bot_input_ids = torch.cat([chat_history_ids, input_ids], dim=-1) if step > 0 else input_ids
        # generate a bot response
        chat_history_ids = model.generate(
            bot_input_ids,
            max_length=1000,
            do_sample=True,
            top_p=0.95,
            top_k=0,
            temperature=0.75,
            pad_token_id=tokenizer.eos_token_id
        )
        #print the output
        output = tokenizer.decode(chat_history_ids[:, bot_input_ids.shape[-1]:][0], skip_special_tokens=True)
        print(f"DialoGPT: {output}")

    # chatting 5 times with nucleus & top-k sampling & tweaking temperature & multiple
    # sentences
    for step in range(arg):
        # take user input
        text = input(">> You:")
        # encode the input and add end of string token
        input_ids = tokenizer.encode(text + tokenizer.eos_token, return_tensors="pt")
        # concatenate new user input with chat history (if there is)
        bot_input_ids = torch.cat([chat_history_ids, input_ids], dim=-1) if step > 0 else input_ids
        # generate a bot response
        chat_history_ids_list = model.generate(
            bot_input_ids,
            max_length=1000,
            do_sample=True,
            top_p=0.95,
            top_k=50,
            temperature=0.75,
            num_return_sequences=5,
            pad_token_id=tokenizer.eos_token_id
        )
        #print the outputs
        for i in range(len(chat_history_ids_list)):
          output = tokenizer.decode(chat_history_ids_list[i][bot_input_ids.shape[-1]:], skip_special_tokens=True)
          print(f"DialoGPT {i}: {output}")
        choice_index = int(input("Choose the response you want for the next input: "))
        chat_history_ids = torch.unsqueeze(chat_history_ids_list[choice_index], dim=0)

    pass

def main(args):
    Testing_chat(5)
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main(sys.argv))
