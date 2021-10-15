import time
import google.api_core.exceptions
from google.cloud import bigquery


def execute_with_retry(client, label, sql, parallelize = False):
    """Run a BigQuery SQL string with a label.

    Three retries with incremental backoff if BiqQuery returns 'retry-able errors'
    (see https://googleapis.dev/python/bigquery/latest/_modules/google/api_core/retry.html),
    any other error is re-raised.

    Parameters
    ----------
    client : bigquery.Client
        with credentials, project, and default_query_job_config already defined
    label : str
        additional label to add to job
    sql : str
        SQL to run
    """
    retry_delay = [30, 60, 90]
    start = time.time()
    while len(retry_delay) >= 0:
        try:
            query = start_query(client, label, sql)
            (results, mb_billed) = get_query_results(query, client, start, retry_delay)
            print(f"COMPLETED ({time.time() - start} seconds, {3 - len(retry_delay)} retries, {mb_billed} MBs) - {label}")

        except (google.api_core.exceptions.InternalServerError,
                google.api_core.exceptions.TooManyRequests,
                google.api_core.exceptions.ServiceUnavailable) as err:
            if len(retry_delay) > 0:
                t = retry_delay.pop(0)
                print(f"Error {err} running query {label}, sleeping for {t}")
                time.sleep(t)
            else:
                raise err
        except Exception:
            raise


def start_query(client, label, sql):
    query_label = label.replace(" ", "-").strip().lower()
    existing_labels = client._default_query_job_config.labels
    job_labels = existing_labels
    job_labels["gvs_query_name"] = query_label
    job_config = bigquery.QueryJobConfig(labels=job_labels)
    query = client.query(sql, job_config=job_config)
    print(f"STARTING - {label} (jobid: {query.job_id})")

    return query


def get_query_results(query, client, start, retry_delay, label):
    results = query.result()
    job = client.get_job(query.job_id)
    mb_billed = int(0 if job.total_bytes_billed is None else job.total_bytes_billed) / (
                1024 * 1024)
    # print(
    #     f"COMPLETED ({time.time() - start} seconds, {3 - len(retry_delay)} retries, {mb_billed} MBs) - {label}")

    return (results, mb_billed)


def utf8len(s):
    return len(s.encode('utf-8'))