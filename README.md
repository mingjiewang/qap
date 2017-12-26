# [![Build Status](https://travis-ci.org/JhuangLab/annovarR.svg)](https://travis-ci.org/JhuangLab/annovarR) [![License](https://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](https://en.wikipedia.org/wiki/MIT_License) [![codecov](https://codecov.io/github/JhuangLab/annovarR/branch/master/graphs/badge.svg)](https://codecov.io/github/JhuangLab/annovarR) 

QAP
==============
[QAP](https://github.com/mingjiewang/qap) (Quasispecies Analysis Package) is an integrated software to analyze viral quasispecies (QS) high through-put sequencing data, including next generation sequencing (NGS) and third generation sequencing (TGS) data. 

Virus community, also known as Quasispecies (QS) is highly related to pathogenesis of viral infectious diseases. Recent development of high through-put sequencing have dramatically lowered the cost and labor of QS detection, yet making computational analysis a major limiting step and an enormous challenge. There is an urgent need for an integrated workflow combining different processing steps in quasispecies studies to discover clinical significance underlying virus populations that could be used on a daily basis by clinicians and virologists. That's why we developed QAP, a powerful all-in-one software to solve the problem.

There are 41 tools included in QAP till now, and all these tools are classified into 6 categories. 
- Raw data preprocessing
- Sequences manipulations
- Quasispecies characterization
- Multiple samples comparison 
- Useful tools
- Visualization and plots

After you download or cloned the source code, use QAP with commands below.
```bash
# Uncompress the tarball file and enter the directory
$ tar xvzf qap.tar.gz
$ cd qap/
# Check dependencies 
$ ./configure
# Install missing dependencies
$ ./autoInstall
# The step of installation and configuration will take some time. Next, add QAP to system PATH.
$ echo "export PATH=`pwd`:$PATH" >> ~/.bash_profile
$ source ~/.bash_profile 
# Now, enjoy! 
$ qap
```

QAP does not only provide command line tools, but also provide a pretty GUI and web applications for users who are uncomfortable with command lines. For more information, please visit [QAP website](http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/).


## Requirements
Referring to hardware requirements for running QAP, there aren't any indispensable ones. If there must be some to be listed, maybe a machine with memory >= 2G, physical threads >= 4 and as much as free disk space, which will give users a stable and faster computing performance. Referring to operating systems, as most of QAP is written Perl, R and Java, which are all cross-platform programming languages. So, theoretically, if you get Perl, R and Java ready on your system, QAP can run on Linux, Windows and Mac OS. However, some of the third-party softwares used by QAP can only run on Linux/Unix based systems, thus making Linux/Unix the recommended system.

We have tested QAP on CentOS, Ubuntu/Debian, RedHat and MacOS, and they all run successfully. Windows OS is not tested yet, but some tools which doesn't depends on Linux based third party tools might run successfully. If you are stucked with Windows OS, please try [Cygwin](http://www.cygwin.com/) or virtual machine softwares, e.g. [VMware](https://www.vmware.com/), which is extraordinary.

For more information about installation and requirement details, please go to [QAP website](http://bioinfo.rjh.com.cn/labs/jhuang/tools/qap/installation/).


## Basic Usage
All tools could be invoked from the wrapped main program. Thus, the usage of QAP is quite easy. Just follow the command below.
```bash
# General command
$ qap \[Tool name\] \[Arguments\]
# Help information
$ qap -h or qap \[Tool name\] -h
# Initiate GUI
$ qap -g

## Docker
We also provide for Docker image for full package of QAP. Download it from [docker hub](https://hub.docker.com/r/mingjiewang/qap/).

```bash
docker pull mingjiewang/qap
docker run -it -v /tmp/db:/tmp/db -v /tmp/input:/tmp/input mingjiewang/qap /usr/bin/bash
```

## How to contribute?

Please fork the [GitHub QAP repository](https://github.com/mingjiewang/qap), modify it, and submit a pull request to us. 

## Maintainer

[Mingjie Wang](https://github.com/mingjiewang/)
Email: [huzai920621@126.com](huzai920621@126.com)
Affiliation: [Ruijin Hospital, Shanghai Jiaotong University, School of Medicine](http://bioinfo.rjh.com.cn/labs/jhuang/)

## License

Main program
[GPL3](https://www.gnu.org/licenses/gpl-3.0.en.html)

R package
[MIT](https://en.wikipedia.org/wiki/MIT_License)

Related Other Resources:
[Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](https://creativecommons.org/licenses/by-nc-nd/4.0/)

