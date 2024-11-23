# Partial Key Exposure Attack on Common Prime RSA

## Introduction

This is a Python implementation of lattice-based attack proposed in **Partial Key Exposure Attack on Common Prime RSA**[^PKEACPRSA].

## Requirements

- [**SageMath**](https://www.sagemath.org/) 9.5 with Python 3.10

You can check your SageMath Python version using the following command:

```commandline
$ sage -python --version
Python 3.10.12
```

Note: If your SageMath Python version is older than 3.9.0, some features in given scripts might not work.

## Usage

The standard way to run the attack with the specific parameters requires passing them as command line arguments `sage -python attack.py <modulus_bit_length> <gamma> <delta> <delta_MSB> <delta_LSB> <t>`. For instance, to run the attack with known MSB, one can run `sage -python attack.py 512 0.48 0.49 0.146 0 6`:

```commandline
PKEA_CPRSA$ sage -python attack.py 512 0.48 0.49 0.146 0 6
The known parameters:
N = 8596626161125351265918067215698641604902111895118306034020399918319803670313027894565152569925510291747556343985793705850786786668155563467345442640857183
e = 57032714234584447497349651007464660310805132876268858972844557181800976153062051
MSB = 10204009548865133291898
LSB = 0
Found primes:
p = 115524498676266928971919668873536269839218903775786868017993964322112703151403
q = 74413879823149764174098202759088832736983665578912343707472352942437960847261
The attack costs 2.317 seconds...
```

For instance, to run the attack with known LSB, one can run `sage -python attack.py 512 0.45 0.327 0 0.041 6`:

```commandline
PKEA_CPRSA$ sage -python attack.py 512 0.45 0.327 0 0.041 6
The known parameters:
N = 10495794222957160392060867587561106048105965233902123090260754390580230944016862332583447513641859880873437847128423345650038838400618092975924447922567415
e = 225885202236748928193001219247116477355773358932234999575234503727193539384443494175
MSB = 0
LSB = 552787
Found primes:
p = 106779923744049723916604309849497348564809958516947241386309790796916609466555
q = 98293704049793677328943859612046838146083172993515983324141174110404824888053
The attack costs 3.119 seconds...
```

For instance, to run the attack with known MSB and LSB, one can run `sage -python attack.py 512 0.44 0.51 0.102 0.145 7`:

```commandline
PKEA_CPRSA$ sage -python attack.py 512 0.44 0.51 0.102 0.145 7
The known parameters:
N = 8839575927486377399313335504884532730183359737276081936109573969770272891353277736452864966051875897027994664843664275681646534012606294195799914064294555
e = 40734180827005069612639483663245354659561607793436239984910141100706735725118784960839
MSB = 4402457607287863
LSB = 1354287992537350484407
Found primes:
p = 101100121868760973824719397742779520541462207179316918148347710918522381564251
q = 87433880039839266172875094503209430896928447393410876733098929640535840698305
The attack costs 19.070 seconds...
```

## Notes

All the details of the numerical attack experiments are recorded in the `attack.log` file.

[^PKEACPRSA]: Zheng M., "Partial Key Exposure Attack on Common Prime RSA"
