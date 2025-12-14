#### ArXiv's PDF Link: [https://arxiv.org/abs/2503.05136](https://arxiv.org/abs/2503.05136)

#### Dynamic Website Version: [https://fhetextbook.github.io](https://fhetextbook.github.io)

-  Please post any bugs or errors regarding the draft to the Issues board or create a pull request.

## Python Demo FHE Library Quickstart

#### Installation
```
pip3 install numpy fastcore sympy
```

### BFV Demo Library
```
usage: bfv_example.py [-h] [--encoding] [--encrypt] [--add-cipher-cipher] 
                    [--add-cipher-plain] [--mult-cipher-cipher] [--mult-cipher-plain] \
                    [--rotate] [--keyswitch] [--bootstrapping] [--random] [--all]

options:
  -h, --help            show this help message and exit
  --encoding            Encoding test
  --random              A bulk of random tests
  --encrypt             Encrytion/decryption test
  --add-cipher-cipher   Cipher-cipher addition test
  --add-cipher-plain    Cipher-plain addition test
  --mult-cipher-cipher  Cipher-cipher multiplication test
  --mult-cipher-plain   Cipher-plain multiplication test
  --rotate              Rotation test
  --keyswitch           Key switch test
  --bootstrapping       Bootstrapping test
  --all                 All test
```
| Test Operations | Command |
|---|---|
| Encode | `python3 bfv_example.py --encoding` |
| Encrypt | `python3 bfv_example.py --encryption` |
| Add Cipher-Cipher | `python3 bfv_example.py --add-cipher-cipher` |
| Add Cipher-Plain | `python3 bfv_example.py --add-cipher-plain` |
| Multiply Cipher-Plain | `python3 bfv_example.py --mult-cipher-plain` |
| Multiply Cipher-Cipher | `python3 bfv_example.py --mult-cipher-cipher` |
| Key Switching | `python3 bfv_example.py --keyswitch` |
| Rotate | `python3 bfv_example.py --rotate` |
| Random 1000 'Encrypt' Tests | `python3 bfv_example.py --random --encrypt` |
| Random 1000 'Rotate' Tests | `python3 bfv_example.py --random --rotate` |
| Random 1000 Tests | `python3 bfv_example.py --all` |