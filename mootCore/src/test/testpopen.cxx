#include <cstdio>
#include <cstdlib>

int main(void) {
  FILE *stream;
  //  unsigned int crc, sz;
  char crc[20], sz[20];
  stream = popen("cksum ../../xml/mood_test.xml", "r");
  // stream = popen("cksum ../../xml/not_a_file.xml", "r");

  if (!stream) {
    printf("whoops\n");
    exit(1);
  }


  //  fscanf(stream, "%i %i", &crc, &sz);
  int nread = fscanf(stream, "%s %s", crc, sz);
  
  if (nread != 2) {
    fprintf(stderr, "Couldn't read from stream \n");
    exit(1);
  }
  if (ferror(stream))  {
    fprintf (stderr, "Output to stream failed.\n");
    exit (1);
  }

  printf("File crc was %s and size was %s\n", crc, sz);

  return 0;
}
