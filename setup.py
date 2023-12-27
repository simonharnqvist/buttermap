from setuptools import setup
setup(
    packages=["buttermap"],
    name='buttermap',
    version='0.0.2',
    install_requires=["ruffus"],
    entry_points={
        'console_scripts': [
            'buttermap=buttermap.buttermap:etode'
        ]
    }
)
