from setuptools import setup
setup(
    packages=["buttermap"],
    name='buttermap',
    version='0.0.1',
    entry_points={
        'console_scripts': [
            'buttermap=buttermap.buttermap:main'
        ]
    }
)