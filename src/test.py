from subprocess import PIPE, Popen

p = Popen(['program', 'arg1'], stdin=PIPE, stdout=PIPE, stderr=PIPE)
output, err = p.communicate(b"input data that is passed to subprocess' stdin")
rc = p.returncode
