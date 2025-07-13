import os, requests, streamlit as st, smtplib
from email.mime.text import MIMEText
from datetime import datetime, timedelta
from typing import List, Dict
from Bio import Entrez
from langchain_openai import AzureChatOpenAI

# ── ENV VARIABLES (set in Streamlit Cloud ▸ Secrets) ───────────
AZURE_ENDPOINT       = os.environ["AZURE_OPENAI_ENDPOINT"]
AZURE_KEY            = os.environ["AZURE_OPENAI_KEY"]

BREVO_SMTP_SERVER    = "smtp-relay.brevo.com"
BREVO_SMTP_PORT      = 587
BREVO_SMTP_USER      = os.environ["BREVO_SMTP_USER"]     # e.g. 91ff83001@smtp-brevo.com
BREVO_SMTP_PASS      = os.environ["BREVO_SMTP_PASS"]     # SMTP password you tested
BREVO_SENDER_EMAIL   = os.environ["BREVO_SENDER_EMAIL"]  # verified sender (e.g. microbiomecollab@gmail.com)
BREVO_SENDER_NAME    = os.environ["BREVO_SENDER_NAME"]   # display name

Entrez.email = "you@example.com"  # PubMed contact email

llm = AzureChatOpenAI(
    azure_endpoint     = AZURE_ENDPOINT,
    api_key            = AZURE_KEY,
    azure_deployment   = "gpt-4o",
    openai_api_version = "2024-12-01-preview",
    temperature        = 0,
)

# ── HELPERS ────────────────────────────────────────────────────
def fetch_pubmed(term: str, max_results: int, lookback: int | None) -> List[Dict]:
    search_args = dict(db="pubmed", term=term, retmax=max_results)
    if lookback:
        since = (datetime.utcnow() - timedelta(days=lookback)).strftime("%Y/%m/%d")
        search_args.update(dict(mindate=since, datetype="pdat"))
    ids = Entrez.read(Entrez.esearch(**search_args))["IdList"]
    if not ids: return []
    docs = Entrez.read(Entrez.esummary(db="pubmed", id=",".join(ids)))
    return [
        {"title": doc["Title"][:250],
         "date":  doc.get("PubDate", "N/A"),
         "url":   f"https://pubmed.ncbi.nlm.nih.gov/{doc['Id']}/"}
        for doc in docs
    ]

def fetch_trials(condition: str, keyword: str, max_results: int) -> List[Dict]:
    r = requests.get("https://clinicaltrials.gov/api/v2/studies", params={
        "query.cond": condition,
        "query.term": keyword,
        "pageSize":   max_results,
        "format":     "json"
    }, timeout=20)
    r.raise_for_status()
    trials = r.json().get("studies", [])
    out = []
    for s in trials:
        ps = s["protocolSection"]
        ident  = ps["identificationModule"]
        status = ps["statusModule"]
        phase  = ps.get("designModule", {}).get("phaseList", {}).get("phase", ["N/A"])[0]
        out.append({
            "title": ident.get("officialTitle") or ident.get("briefTitle"),
            "nct":   ident["nctId"],
            "phase": phase,
            "start": status.get("startDateStruct", {}).get("startDate", "N/A"),
            "url":   f"https://clinicaltrials.gov/study/{ident['nctId']}",
        })
    return out

def bullet(text: str) -> str:
    prompt = ("Summarize in ≤30 words for oncologists:\n\n" + text + "\n\nBullet:")
    return llm.invoke(prompt).content.strip()

def send_email(to_email: str, subject: str, body: str):
    msg = MIMEText(body)
    msg["Subject"] = subject
    msg["From"] = f"{BREVO_SENDER_NAME} <{BREVO_SENDER_EMAIL}>"
    msg["To"] = to_email
    with smtplib.SMTP(BREVO_SMTP_SERVER, BREVO_SMTP_PORT) as srv:
        srv.starttls()
        srv.login(BREVO_SMTP_USER, BREVO_SMTP_PASS)
        srv.send_message(msg)

# ── STREAMLIT UI ───────────────────────────────────────────────
st.title("🧬 Oncology Digest Generator")

topic          = st.text_input("Search term", "ctDNA bladder cancer")
recipient      = st.text_input("Recipient email")
max_papers     = st.slider("Max papers", 1, 100, 8)
max_trials     = st.slider("Max trials", 1, 50, 5)
lookback_opt   = st.selectbox("PubMed window", ("7 days", "30 days", "All time"))
lookback_days  = {"7 days": 7, "30 days": 30, "All time": None}[lookback_opt]

if st.button("Generate & Send"):
    with st.spinner("Collecting literature…"):
        condition = "bladder cancer"
        keyword   = "ctDNA" if "ctDNA" in topic else topic
        papers = fetch_pubmed(topic, max_papers, lookback_days)
        trials = fetch_trials(condition, keyword, max_trials)

    for p in papers:  p["summary"] = bullet(p["title"])
    for t in trials:  t["summary"] = bullet(t["title"])

    today = datetime.now().strftime("%d %b %Y")
    email_lines = [
        f"Subject: Oncology Digest – {topic} ({today})",
        "",
        "Dear Dr. Bukavina,",
        "",
        f"Digest on *{topic}*: ",
        "",
        "📚 Papers:"
    ] + (
        sum([[f"• {p['summary']}", f"  ({p['date']}) {p['url']}"] for p in papers], [])
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
    digest = "\n".join(email_lines)

    st.text_area("Digest preview", digest, height=400)

    if recipient:
        try:
            send_email(recipient, f"Oncology Digest – {topic}", digest)
            st.success(f"Email sent to {recipient}")
        except Exception as e:
            st.error(f"Email failed: {e}")
    else:
        st.warning("Enter a recipient email to send.")
