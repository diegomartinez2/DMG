#!/bin/bash
#transformar de video mp4 a audio mp3 con ffmpeg
ffmpeg -i Musica.mp4 -f mp3 -ab 192000 -vn Musica.mp3
