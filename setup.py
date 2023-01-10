from setuptools import setup

setup(name='biotools_dxy',
      version='0.1',
      description='self-defined functions for bioinformatics',
      long_description="self-defined functions for bioinformatics",
      classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3'
        ],
      keywords='bioinformatics tools',
      url='https://github.com/yang-dongxu/biotools_dxy',
      author='Dongxu Yang',
      author_email='yang_dongxu@qq.com',
      license='MIT',
      packages=['biotools_dxy'],
      install_requires=[
        'numpy',
        'pandas',
        'bioframe',
        "pyranges"
      ],
      include_package_data=True,
      zip_safe=False)
