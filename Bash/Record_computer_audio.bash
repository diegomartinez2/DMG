#Not tested yet
arecord -f cd -D plughw:1,0 | ffmpeg -i - -vn -ar 44100 -ac 2 -b:a 192k output.mp3
#this should record the audio output from a linux pc computer into a 'output.mp3' file.
