# Build PDF locally

See [docs][docs] for detailed instructions.

- Install docker
- Navigate to directory which contains paper.md
- Execute the command:

```
sudo docker run --rm \
    --volume $PWD:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/paperdraft
```

[docs]: https://joss.readthedocs.io/en/latest/submitting.html#example-paper-and-bibliography
