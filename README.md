# pspTSGen

 Currently not working.

### Installation & Docker

Install the dependencies and start the server.

```sh
$ git clone https://github.com/rackerm4/pspTSGen

$ docker build -t pspTSGen .

$ docker run pspTSGen --profile default --num_runs <N> --schema newick --output data

```

### Arguments


| Arg | ... |
| ------ | ------ |
| --profile | |
| --num_runs | |
| --schema | |
| --output | |

### Requirements

* DendroPy==4.4.0
* future==0.18.2
* PyYAML==5.3.1

