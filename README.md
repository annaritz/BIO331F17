# BIO331F17
Group project to identify potential regulators of NMII and Fog pathway members in a fly interactome.  This README is written using [MarkDown](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet).

## Quick Start

This is a quick start for the toy script.  Run the command

```
python3 HelloWorld.py
```

## Common GitHub Commands

- `git pull`: "pull" the version on GitHub to your local machine. **Always do this before you start working.**
- `git status`: print the status of all the tracked and untracked files.
- `git add <file>`: add a file ("stage" it) for a commit. **Do this occaisonally when you are working.**
- `git add -u`: stage **all** tracked files for a commit. **Do this occaisonally when you are working.**
- `git commit -m "<message>"`: commit the staged files on your local machine. **Do this occaisonally when you are working.**
- `git push origin head`: "push" the commits to GitHub. **Always do this when you're done working.**

## GitHub Tutorial

1. You should already have a GitHub account and have been added as collaborators to this repository.  Ensure that you have git installed by typing the following in a terminal:

```
git --version
```

If you do not have `git` installed, [install it according to your OS](https://git-scm.com/downloads).

2. "Clone" the project onto your local machine.  Change directories to where you want to esablish the project.  Click the green "Clone or Download" button in the upper right corner and copy the link.  You can then paste the link when typing this command:

```
git clone https://github.com/annaritz/BIO331F17.git
```

You can then change directories to the newly-established `BIO331F17`.

3.  Open `HelloWorld.py` and confirm that it looks the same as the version on GitHub.  

4. Make a change to `HelloWorld.py`.  Type

```
git status
```

and see what has changed.

5. Now, you will **commit** this change.  Commits allow you to make incremental changes to a program and be able to "roll back" to a previous commit at any time.  There are two types of commits: (a) committing the change your local machine and (b) **pushing** those commits to the repository on GitHub.  

5(a).  Commit the change to your local machine.  First, you must **stage** the files for commit by telling git what to keep track of.  This is useful when you are modifying a bunch of files but only want to commit a subset of them.

```
git add HelloWorld.py
```

Type `git status` and observe the change.  Then, you can commit the change:

```
git commit -m "making small change in HelloWorld"
```

The `-m` argument is the *message*, which should be descriptive.  If you forget the `-m` argument, you may get an error or be redirected to a text editor.

5(b). Push the commit to GitHub:

```
git push origin head
```

This is saying that we will push the commits from the local machine (the **origin**) to the repository on GitHub (the **head**).

6. You might have gotten an error at this point!  Someone else might have pushed their commits before you did. In that case, you need to **pull** their changes to your local machine.

```
git pull
```

You can then try pushing your commits.


