# Contributing

You are welcome to:

* submit suggestions and bug-reports at: <https://github.com/sgibb/topdown/issues>
* send a pull request on: <https://github.com/sgibb/topdown/compare>
* compose an e-mail to: <mail@sebastiangibb.de>

# How to contribute code

[Fork](https://help.github.com/articles/fork-a-repo/), then clone the repository:

    git clone git@github.com:your-username/topdown.git

Make sure all checks and tests pass:

    R CMD build topdown && CMD check --as-cran --no-stop-on-test-error topdown_*.tar.gz

Make your changes and write tests ...

Ensure all checks and tests pass:

    R CMD build topdown && CMD check --as-cran --no-stop-on-test-error topdown_*.tar.gz

Push to your fork and submit a [pull request](https://help.github.com/articles/about-pull-requests/)

Waiting for our response. We may suggest some changes and/or improvements.

If you want to increase the chance that your pull request is accepted:

* Use the
  [Bioconductor Coding Style](https://www.bioconductor.org/developers/how-to/coding-style/)
  with the exception that we prefer 2 spaces instead of 4 for indention.
* Write [good commit messages](https://robots.thoughtbot.com/5-useful-tips-for-a-better-commit-message).
* Write unit tests using the `testthat` package.
