all:
	gcc main.c tga.c model.c -o main -lm -D DEBUG
clean:
	rm main
