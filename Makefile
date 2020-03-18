CC = gcc
RM = rm -f
CFLAGS = -I.

SRCS = $(wildcard *.c)
OBJS = $(SRCS:.c=.o)
TARGET_EXE = fours.exe

.PHONY: all
all: $(TARGET_EXE)

$(TARGET_EXE): $(OBJS)
	$(CC) -o $@ $^ -lgmp

.PHONY: clean
clean:
	-$(RM) ${TARGET_EXE} ${OBJS} $(wildcard *.*~)

.PHONY: protobuf
protobuf:
	cd res; protoc --c_out=.. CredentialSchema.proto ProofHashData.proto
