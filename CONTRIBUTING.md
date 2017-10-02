# Contributing

You are welcome to:

* submit suggestions and bug-reports at: <https://github.com/sgibb/topdownr/issues>
* send a pull request on: <https://github.com/sgibb/topdownr/compare>
* compose an e-mail to: <mail@sebastiangibb.de>

# How to contribute code

[Fork](https://help.github.com/articles/fork-a-repo/), then clone the repository:

    git clone git@github.com:your-username/topdownr.git

Make sure all checks and tests pass:

    R CMD build topdownr && CMD check --as-cran --no-stop-on-test-error topdownr_*.tar.gz

Make your changes and write tests ...

Ensure all checks and tests pass:

    R CMD build topdownr && CMD check --as-cran --no-stop-on-test-error topdownr_*.tar.gz

Push to your fork and submit a [pull request](https://help.github.com/articles/about-pull-requests/)

Waiting for our response. We may suggest some changes and/or improvements.

If you want to increase the chance that your pull request is accepted:

* Use the
  [Bioconductor Coding Style](https://www.bioconductor.org/developers/how-to/coding-style/).
* Write [good commit messages](https://robots.thoughtbot.com/5-useful-tips-for-a-better-commit-message).
* Write unit tests using the `testthat` package.
