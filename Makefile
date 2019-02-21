LFLAGS = -lm
all:
	gcc main.c tga.c model.c -o main
clean:
	rm main