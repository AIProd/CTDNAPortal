import os, re, requests, streamlit as st, smtplib
from email.mime.text import MIMEText
from datetime import datetime, timedelta
from typing import List, Dict
from Bio import Entrez
from langchain_openai import AzureChatOpenAI

# ── SECRETS (Streamlit Cloud ▸ Settings ▸ Secrets) ────────────
AZURE_ENDPOINT       = os.environ["AZURE_OPENAI_ENDPOINT"]
AZURE_KEY            = os.environ["AZURE_OPENAI_KEY"]
BREVO_SMTP_USER      = os.environ["BREVO_SMTP_USER"]
BREVO_SMTP_PASS      = os.environ["BREVO_SMTP_PASS"]
BREVO_SENDER_EMAIL   = os.environ["BREVO_SENDER_EMAIL"]
BREVO_SENDER_NAME    = os.environ["BREVO_SENDER_NAME"]

Entrez.email = "you@example.com"  # PubMed contact

llm = AzureChatOpenAI(
    azure_endpoint     = AZURE_ENDPOINT,
    api_key            = AZURE_KEY,
    azure_deployment   = "gpt-4o",
    openai_api_version = "2024-12-01-preview",
    temperature        = 0,
)

# ── Helper functions ───────────────────────────────────────────
def is_valid_email(addr: str) -> bool:
    return bool(re.match(r"[^@]+@[^@]+\.[^@]+", addr))

def fetch_pubmed(term: str, n: int, lookback: int | None) -> List[Dict]:
    args = dict(db="pubmed", term=term, retmax=n)
    if lookback:
        since = (datetime.utcnow() - timedelta(days=lookback)).strftime("%Y/%m/%d")
        args.update(mindate=since, datetype="pdat")
    ids = Entrez.read(Entrez.esearch(**args))["IdList"]
    if not ids:
        return []
    docs = Entrez.read(Entrez.esummary(db="pubmed", id=",".join(ids)))
    return [
        {
            "title": d["Title"][:250],
            "date":  d.get("PubDate", "N/A"),
            "url":   f"https://pubmed.ncbi.nlm.nih.gov/{d['Id']}/",
        }
        for d in docs
    ]

def fetch_trials(cond: str, kw: str, n: int) -> List[Dict]:
    r = requests.get("https://clinicaltrials.gov/api/v2/studies", params={
        "query.cond": cond,
        "query.term": kw,
        "pageSize":   n,
        "format":     "json"},
        timeout=20)
    r.raise_for_status()
    out = []
    for s in r.json().get("studies", []):
        ps   = s["protocolSection"]
        ident = ps["identificationModule"]
        phase = ps.get("designModule", {}).get("phaseList", {}).get("phase", ["N/A"])[0]
        out.append({
            "title": ident.get("officialTitle") or ident["briefTitle"],
            "nct":   ident["nctId"],
            "phase": phase,
            "start": ps["statusModule"].get("startDateStruct", {}).get("startDate", "N/A"),
            "url":   f"https://clinicaltrials.gov/study/{ident['nctId']}"
        })
    return out

def bullet(txt: str) -> str:
    prompt = ("Summarize in ≤30 words for oncologists:\n\n" + txt + "\n\nBullet:")
    return llm.invoke(prompt).content.strip().lstrip("- ").strip()

def send_email(to_addr: str, subj: str, body: str):
    msg = MIMEText(body)
    msg["Subject"] = subj
    msg["From"]    = f"{BREVO_SENDER_NAME} <{BREVO_SENDER_EMAIL}>"
    msg["To"]      = to_addr
    with smtplib.SMTP("smtp-relay.brevo.com", 587) as srv:
        srv.starttls()
        srv.login(BREVO_SMTP_USER, BREVO_SMTP_PASS)
        srv.send_message(msg)

# ── Streamlit UI ───────────────────────────────────────────────
st.title("🧬 Oncology Digest Generator")

topic     = st.text_input("Search term", "ctDNA bladder cancer")
name      = st.text_input("Recipient name", placeholder="Enter recipient name")
email     = st.text_input("Recipient email", placeholder="Enter a valid email")
papers_n  = st.slider("Max papers", 1, 100, 8)
trials_n  = st.slider("Max trials", 1, 50, 5)
lookback  = st.selectbox("PubMed window", ("7 days", "30 days", "All time"))
days_map  = {"7 days": 7, "30 days": 30, "All time": None}
days      = days_map[lookback]

if st.button("Generate & Send"):

    # ── Validation ─────────────────────────────────────────────
    if not topic.strip():
        st.error("Please enter a search term.")
        st.stop()
    if not name.strip():
        st.error("Please enter the recipient's name.")
        st.stop()
    if not is_valid_email(email):
        st.error("Please enter a valid email address.")
        st.stop()

    progress = st.progress(0)
    with st.status("Generating digest…", expanded=True) as status:
        status.write("🔍 Fetching PubMed articles…")
        papers = fetch_pubmed(topic, papers_n, days); progress.progress(20)

        status.write("🔍 Fetching clinical trials…")
        trials = fetch_trials("bladder cancer", topic, trials_n); progress.progress(40)

        status.write("🧠 Summarizing with GPT…")
        for p in papers:
            p["summary"] = bullet(p["title"])
        progress.progress(60)
        for t in trials:
            t["summary"] = bullet(t["title"])
        progress.progress(80)

        today = datetime.now().strftime("%d %b %Y")
        body_lines = [
            f"Dear {name},",
            "",
            f'Here is your digest on "{topic}":',
            "",
            "📚 Papers:"
        ] + (
            sum([[f"• {p['summary']}",
                  f"  ({p['date']}) {p['url']}"] for p in papers], [])
            or ["  – None found."]
        ) + [
            "",
            "🧪 Trials:"
        ] + (
            sum([[f"• {t['summary']} [{t['phase']}]",
                  f"  NCT {t['nct']} — Start {t['start']}",
                  f"  {t['url']}"] for t in trials], [])
            or ["  – None found."]
        ) + ["", "Regards,", "Oncology AI"]

        digest_body = "\n".join(body_lines)
        subject = f"Oncology Digest – {topic} ({today})"
        progress.progress(90)

        try:
            send_email(email, subject, digest_body)
            status.success("✅ Email sent!")
        except Exception as exc:
            status.error(f"Email failed: {exc}")

        progress.progress(100)
        status.write("🚀 Done!")

    st.text_area("Digest preview", digest_body, height=400)
