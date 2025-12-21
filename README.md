#### ArXiv's PDF Link: [https://arxiv.org/abs/2503.05136](https://arxiv.org/abs/2503.05136)

#### Dynamic Website Version: [https://fhetextbook.github.io](https://fhetextbook.github.io)

-  Please post any bugs or errors regarding the draft to the Issues board or create a pull request.

## Python Demo FHE Library Quickstart


This library is provided for an educational purpose.



### Installation




```
cd 'source code'
pip3 install numpy fastcore sympy
```


### TFHE Library
```
usage: tfhe_example.py [-h] [--encrypt] [--add-cipher-cipher] 
                    [--add-cipher-plain] [--mult-cipher-cipher] [--mult-cipher-plain] \
                    [--keyswitch] [--random] [--all]

options:
  -h, --help            show this help message and exit
  --random              A bulk of random tests
  --encrypt             Encrytion/decryption test
  --add-cipher-cipher   Cipher-cipher addition test
  --add-cipher-plain    Cipher-plain addition test
  --mult-cipher-cipher  Cipher-cipher multiplication test
  --mult-cipher-plain   Cipher-plain multiplication test
  --keyswitch           Key switch test
  --all                 All test
```
| Test Operations | Command |
|---|---|
| Encrypt | `python3 tfhe_example.py --encryption` |
| Add Cipher-Cipher | `python3 tfhe_example.py --add-cipher-cipher` |
| Add Cipher-Plain | `python3 tfhe_example.py --add-cipher-plain` |
| Multiply Cipher-Plain | `python3 tfhe_example.py --mult-cipher-plain` |
| Multiply Cipher-Cipher | `python3 tfhe_example.py --mult-cipher-cipher` |
| Key Switching | `python3 tfhe_example.py --keyswitch` |
| Random 1000 'Encrypt' Tests | `python3 tfhe_example.py --random --encrypt` |
| Random 1000 'Rotate' Tests | `python3 tfhe_example.py --random --rotate` |
| Random 1000 Any Tests | `python3 tfhe_example.py --all` |


### BFV Library
```
usage: bfv_example.py [-h] [--encode] [--encrypt] [--add-cipher-cipher] 
                    [--add-cipher-plain] [--mult-cipher-cipher] [--mult-cipher-plain] \
                    [--rotate] [--keyswitch] [--random] [--all]

options:
  -h, --help            show this help message and exit
  --encode              Encoding test
  --random              A bulk of random tests
  --encrypt             Encrytion/decryption test
  --add-cipher-cipher   Cipher-cipher addition test
  --add-cipher-plain    Cipher-plain addition test
  --mult-cipher-cipher  Cipher-cipher multiplication test
  --mult-cipher-plain   Cipher-plain multiplication test
  --rotate              Rotation test
  --keyswitch           Key switch test
  --all                 All test
```
| Test Operations | Command |
|---|---|
| Encode | `python3 bfv_example.py --encode` |
| Encrypt | `python3 bfv_example.py --encryption` |
| Add Cipher-Cipher | `python3 bfv_example.py --add-cipher-cipher` |
| Add Cipher-Plain | `python3 bfv_example.py --add-cipher-plain` |
| Multiply Cipher-Plain | `python3 bfv_example.py --mult-cipher-plain` |
| Multiply Cipher-Cipher | `python3 bfv_example.py --mult-cipher-cipher` |
| Key Switching | `python3 bfv_example.py --keyswitch` |
| Rotate | `python3 bfv_example.py --rotate` |
| Conjugate | `python3 bfv_example.py --conjugate` |
| Random 1000 'Encrypt' Tests | `python3 bfv_example.py --random --encrypt` |
| Random 1000 'Rotate' Tests | `python3 bfv_example.py --random --rotate` |
| Random 1000 Any Tests | `python3 bfv_example.py --all` |


### CKKS Library
```
usage: ckks_example.py [-h] [--encode] [--encrypt] [--add-cipher-cipher] 
                    [--add-cipher-plain] [--mult-cipher-cipher] [--mult-cipher-plain] \
                    [--rotate] [--conjugate] [--keyswitch] [--bootstrapping] [--random] [--all]

options:
  -h, --help            show this help message and exit
  --encode              Encoding test
  --random              A bulk of random tests
  --encrypt             Encrytion/decryption test
  --add-cipher-cipher   Cipher-cipher addition test
  --add-cipher-plain    Cipher-plain addition test
  --mult-cipher-cipher  Cipher-cipher multiplication test
  --mult-cipher-plain   Cipher-plain multiplication test
  --rotate              Rotation test
  --conjugate           Conjugation test
  --keyswitch           Key switch test
  --all                 All test
```
| Test Operations | Command |
|---|---|
| Encode | `python3 ckks_example.py --encode` |
| Encrypt | `python3 ckks_example.py --encryption` |
| Add Cipher-Cipher | `python3 ckks_example.py --add-cipher-cipher` |
| Add Cipher-Plain | `python3 ckks_example.py --add-cipher-plain` |
| Multiply Cipher-Plain | `python3 ckks_example.py --mult-cipher-plain` |
| Multiply Cipher-Cipher | `python3 ckks_example.py --mult-cipher-cipher` |
| Key Switching | `python3 ckks_example.py --keyswitch` |
| Rotate | `python3 ckks_example.py --rotate` |
| Conjugate | `python3 ckks_example.py --conjugate` |
| Random 1000 'Encrypt' Tests | `python3 ckks_example.py --random --encrypt` |
| Random 1000 'Rotate' Tests | `python3 ckks_example.py --random --rotate` |
| Random 1000 Any Tests | `python3 ckks_example.py --all` |



### BGV Library
```
usage: ckks_example.py [-h] [--encode] [--encrypt] [--add-cipher-cipher] 
                    [--add-cipher-plain] [--mult-cipher-cipher] [--mult-cipher-plain] \
                    [--rotate] [--keyswitch] [--bootstrapping] [--random] [--all]

options:
  -h, --help            show this help message and exit
  --encode              Encoding test
  --random              A bulk of random tests
  --encrypt             Encrytion/decryption test
  --add-cipher-cipher   Cipher-cipher addition test
  --add-cipher-plain    Cipher-plain addition test
  --mult-cipher-cipher  Cipher-cipher multiplication test
  --mult-cipher-plain   Cipher-plain multiplication test
  --rotate              Rotation test
  --keyswitch           Key switch test
  --all                 All test
```
| Test Operations | Command |
|---|---|
| Encode | `python3 bgv_example.py --encoding` |
| Encrypt | `python3 bgv_example.py --encryption` |
| Add Cipher-Cipher | `python3 bgv_example.py --add-cipher-cipher` |
| Add Cipher-Plain | `python3 bgv_example.py --add-cipher-plain` |
| Multiply Cipher-Plain | `python3 bgv_example.py --mult-cipher-plain` |
| Multiply Cipher-Cipher | `python3 bgv_example.py --mult-cipher-cipher` |
| Key Switching | `python3 bgv_example.py --keyswitch` |
| Rotate | `python3 bgv_example.py --rotate` |
| Random 1000 'Encrypt' Tests | `python3 bgv_example.py --random --encrypt` |
| Random 1000 'Rotate' Tests | `python3 bgv_example.py --random --rotate` |
| Random 1000 Any Tests | `python3 bgv_example.py --all` |

## Community Implementations
- For interested readers, a BFV-RNS implementation in SystemVerilog is available at [https://github.com/alibillalhammoud/FHEforEEs](https://github.com/alibillalhammoud/FHEforEEs).
