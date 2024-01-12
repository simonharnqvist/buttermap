from setuptools import setup
setup(
    packages=["buttermap"],
    name='buttermap',
    version='0.0.3',
    install_requires=["ruffus"],
    entry_points={
        'console_scripts': [
            'buttermap=buttermap.buttermap:main'
        ]
    }
)
