# spcount

This package is used to map short reads to multiple genomes with perfect match and report only once for each category.

# Prerequisites

Install bowtie

```
BOWTIE_VERSION="1.2.3"
cd ~; \
  wget https://github.com/BenLangmead/bowtie/releases/download/v${BOWTIE_VERSION}/bowtie-${BOWTIE_VERSION}-linux-x86_64.zip; \
  unzip bowtie-${BOWTIE_VERSION}-linux-x86_64.zip; \
  rm bowtie-${BOWTIE_VERSION}-linux-x86_64.zip
export PATH=$PATH:~/bowtie-${BOWTIE_VERSION}-linux-x86_64
```

# Installation

Install python main package

```
pip install spcount
```

Or you can install from github directly

```
pip install git+git://github.com/shengqh/spcount.git
```

