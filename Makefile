all:
	gcc main.c tga.c model.c -o main -lm
clean:
	rm main
